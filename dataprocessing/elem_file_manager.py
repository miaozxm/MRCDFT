"""
MRCDFT 输出文件 (.elem) 管理器
==============================
统一管理 Proj_*_kern.*.elem 和 Proj_*_R2_2b.*.elem 文件的：
  1. 文件名解析（正则提取元数据）
  2. 文件内容解析（定宽 / 空格分隔）
  3. beta2 交换复制（对称化 kernel 矩阵）
  4. 批量处理 → DataFrame 输出
  5. 增量处理（跳过已处理文件）

用法:
    from elem_file_manager import ElemFileManager

    mgr = ElemFileManager("/home/xizhang/MRCDFT-master/Labs/40Ca/40Ca_10_2026-07-14-09-44-24")
    df = mgr.process_all()          # 处理全部 → DataFrame
    df_kern = mgr.process("kern")   # 只处理 kernel 文件
    df_r2   = mgr.process("R2_2b")  # 只处理 R2 两体文件
"""

import re
import os
import shutil
from dataclasses import dataclass, field
from typing import Optional
import pandas as pd


# ──────────────────────────────────────────────
# 1. 文件类型注册表：pattern + format 定义集中管理
# ──────────────────────────────────────────────

FILE_TYPES = {
    "kern": {
        "pattern": re.compile(
            r'^(?P<prefix>Proj_)'
            r'(?P<m>\d{2})'
            r'(?P<name>[A-Z][a-z]?)'
            r'(?P<filename>_kern)\.'
            r'(?P<dim>\dD)'
            r'(?P<s>_eMax)'
            r'(?P<n0f>\d{2})\.'
            r'(?P<nphi>\d{2})\.'
            r'(?P<nbeta>\d{2})'
            r'(?P<beta2_l>[+-]\d{3})'
            r'(?P<beta3_l>[+-]\d{3})_'
            r'(?P<beta2_r>[+-]\d{3})'
            r'(?P<beta3_r>[+-]\d{3})\.elem'
        ),
        "format": {
            1: {'method': 'space', 'key_point': ['J', 'K1', 'K2', 'Parity', 'N^2/N_KK', 'Z^2/N_KK', 'J^2/N_KK']},
            2: {'method': 'width', 'key_point': [("N_KKr", 15), ('N_KKi', 15), ('H_KKr', 15), ('H_KKi', 15)]},
            3: {'method': 'width', 'key_point': [("Nr", 15), ('Ni', 15), ('Zr', 15), ('Zi', 15)]},
            4: {'method': 'width', 'key_point': [("Q2_KK_12pr", 15), ('Q2_KK_12pi', 15), ('Q_KK_21pr', 15), ('Q_KK_21pi', 15)]},
            5: {'method': 'width', 'key_point': [("E0_KKpr", 15), ('E0_KKpi', 15), ('E0_KKpr', 15), ('E0_KKpi', 15)]},
            6: {'method': 'width', 'key_point': [("Q2_KK_12nr", 15), ('Q2_KK_12ni', 15), ('Q_KK_21nr', 15), ('Q_KK_21ni', 15)]},
            7: {'method': 'width', 'key_point': [("E0_KKnr", 15), ('E0_KKni', 15), ('E0_KKnr', 15), ('E0_KKni', 15)]},
        },
        "lines_per_record": 7,
        "name_fields": ['n0f', 'beta2_l', 'beta3_l', 'beta2_r', 'beta3_r', 'nbeta'],
    },
    "R2_2b": {
        "pattern": re.compile(
            r'^(?P<prefix>Proj_)'
            r'(?P<m>\d{2})'
            r'(?P<name>[A-Z][a-z]?)'
            r'(?P<filename>_R2_2b)\.'
            r'(?P<dim>\dD)'
            r'(?P<s>_eMax)'
            r'(?P<n0f>\d{2})\.'
            r'(?P<nphi>\d{2})\.'
            r'(?P<nbeta>\d{2})'
            r'(?P<beta2_l>[+-]\d{3})'
            r'(?P<beta3_l>[+-]\d{3})_'
            r'(?P<beta2_r>[+-]\d{3})'
            r'(?P<beta3_r>[+-]\d{3})\.elem'
        ),
        "format": {
            1: {'method': 'width', 'key_point': [("total_r2_1B", 15), ('total_r2_2B', 15), ('r2_1&2B', 15)]},
            2: {'method': 'width', 'key_point': [("total_direct", 15), ('total_exchange', 15), ('total_kappa', 15)]},
            3: {'method': 'width', 'key_point': [("o2", 15), ('o5', 15), ('1&2bn', 15)]},
            4: {'method': 'width', 'key_point': [("2bndirect", 15), ('2bnexchange', 15), ('2bnkappa', 15)]},
            5: {'method': 'width', 'key_point': [("o1", 15), ('o3', 15), ('1&2bp', 15)]},
            6: {'method': 'width', 'key_point': [("2bpdirect", 15), ('2bpexchange', 15), ('2bpkappa', 15)]},
            7: {'method': 'width', 'key_point': [("1bpn", 15), ('o4', 15), ('2bpn', 15)]},
            8: {'method': 'width', 'key_point': [("2bpndirect", 15), ('exchange', 15), ('kappa', 15)]},
        },
        "lines_per_record": 8,
        "name_fields": ['n0f', 'beta2_l', 'beta3_l', 'beta2_r', 'beta3_r', 'nbeta'],
    },
}


# ──────────────────────────────────────────────
# 2. 单文件解析器（纯函数，无副作用）
# ──────────────────────────────────────────────

def _char_to_float(char: str):
    """字符串转 float，失败则保留原值（如 Parity 的 '+'）"""
    try:
        return float(char)
    except ValueError:
        return char


def _parse_line(line: str, fmt: dict) -> dict:
    """根据格式描述解析单行，返回 {字段名: 值}"""
    result = {}
    if fmt['method'] == 'space':
        fields = line.split()
        for i, field in enumerate(fields):
            if i < len(fmt['key_point']):
                result[fmt['key_point'][i]] = _char_to_float(field)
    elif fmt['method'] == 'width':
        pos = 0
        for name, width in fmt['key_point']:
            result[name] = _char_to_float(line[pos:pos + width])
            pos += width
    return result


def parse_elem_file(filepath: str, fmt: dict, lines_per_record: int) -> dict:
    """
    解析单个 .elem 文件，返回包含所有字段的总字典。
    假设文件恰好包含 lines_per_record 行有效数据。
    """
    data = {}
    with open(filepath, 'r') as f:
        contents = f.readlines()

    for i, line in enumerate(contents):
        stripped = line.strip()
        if not stripped:
            continue
        line_idx = (i % lines_per_record) + 1
        if line_idx in fmt:
            data.update(_parse_line(line, fmt[line_idx]))

    return data


# ──────────────────────────────────────────────
# 3. 文件名工具（修正了 swap 的位置 bug）
# ──────────────────────────────────────────────

def extract_metadata(filename: str, pattern: re.Pattern) -> Optional[dict]:
    """从文件名提取元数据，beta 值乘以 0.01 转为物理值"""
    match = pattern.match(filename)
    if not match:
        return None

    raw = match.groupdict()
    meta = {}
    for key, val in raw.items():
        if key in ('beta2_l', 'beta3_l', 'beta2_r', 'beta3_r'):
            meta[key] = float(val) * 0.01
        else:
            try:
                meta[key] = float(val)
            except ValueError:
                meta[key] = val
    return meta


def swap_beta2_in_filename(filename: str, pattern: re.Pattern) -> Optional[str]:
    """
    交换文件名中的 beta2_l 和 beta2_r，返回新文件名。
    使用 match.span() 精确定位替换，避免值冲突。
    """
    match = pattern.match(filename)
    if not match:
        return None

    span_l = match.span('beta2_l')  # (start, end) 在原字符串中的位置
    span_r = match.span('beta2_r')
    val_l = match.group('beta2_l')
    val_r = match.group('beta2_r')

    # 按位置从左到右替换（先处理靠前的那个）
    if span_l[0] < span_r[0]:
        new = filename[:span_l[0]] + val_r + filename[span_l[1]:span_r[0]] + val_l + filename[span_r[1]:]
    else:
        new = filename[:span_r[0]] + val_l + filename[span_r[1]:span_l[0]] + val_r + filename[span_l[1]:]
    return new


# ──────────────────────────────────────────────
# 4. 主管理器类
# ──────────────────────────────────────────────

@dataclass
class ElemFileManager:
    """
    .elem 文件管理器：扫描目录、解析内容、交换 beta2、输出 DataFrame。

    Parameters
    ----------
    dirpath : str
        包含 .elem 文件的目录路径
    do_swap_copy : bool
        是否自动复制交换 beta2 后的文件（默认 True）
    skip_existing : bool
        交换复制时，目标文件已存在则跳过（默认 True）
    """
    dirpath: str
    do_swap_copy: bool = True
    skip_existing: bool = True
    _dir_listing: set = field(default_factory=set, init=False)

    def __post_init__(self):
        self._refresh_listing()

    def _refresh_listing(self):
        """刷新目录文件列表缓存"""
        self._dir_listing = set(os.listdir(self.dirpath))

    # ── 单类型处理 ──

    def process(self, file_type: str) -> pd.DataFrame:
        """
        处理指定类型的所有文件，返回 DataFrame。

        Parameters
        ----------
        file_type : str
            "kern" 或 "R2_2b"

        Returns
        -------
        pd.DataFrame
            每行一个文件，列 = 文件名元数据 + 文件内容字段
        """
        if file_type not in FILE_TYPES:
            raise ValueError(f"未知文件类型: {file_type}，可选: {list(FILE_TYPES.keys())}")

        spec = FILE_TYPES[file_type]
        pattern = spec["pattern"]
        fmt = spec["format"]
        lines = spec["lines_per_record"]
        name_fields = spec["name_fields"]

        records = []

        # 先收集所有匹配文件，避免在循环中修改目录
        matched_files = []
        for fname in sorted(self._dir_listing):
            if pattern.match(fname):
                matched_files.append(fname)

        # 批量做交换复制（如果开启）
        if self.do_swap_copy:
            self._batch_swap_copy(matched_files, pattern)

        # 解析每个文件
        for fname in matched_files:
            filepath = os.path.join(self.dirpath, fname)
            try:
                content_data = parse_elem_file(filepath, fmt, lines)
            except Exception as e:
                print(f"[WARN] 解析失败 {fname}: {e}")
                continue

            meta = extract_metadata(fname, pattern)
            # 只保留 name_fields 指定的元数据列（加上 filename 供溯源）
            meta_filtered = {k: meta.get(k) for k in name_fields}
            meta_filtered["filename"] = fname

            records.append({**meta_filtered, **content_data})

        df = pd.DataFrame(records)
        return df

    def _batch_swap_copy(self, filenames: list, pattern: re.Pattern):
        """批量执行 beta2 交换复制"""
        for fname in filenames:
            new_name = swap_beta2_in_filename(fname, pattern)
            if new_name is None or new_name == fname:
                continue
            if self.skip_existing and new_name in self._dir_listing:
                continue
            src = os.path.join(self.dirpath, fname)
            dst = os.path.join(self.dirpath, new_name)
            shutil.copy2(src, dst)
            self._dir_listing.add(new_name)

    # ── 全部处理 ──

    def process_all(self) -> dict[str, pd.DataFrame]:
        """
        处理所有已知文件类型。

        Returns
        -------
        dict
            {"kern": DataFrame, "R2_2b": DataFrame}
        """
        results = {}
        for ftype in FILE_TYPES:
            df = self.process(ftype)
            results[ftype] = df
            print(f"[{ftype}] 处理完成: {len(df)} 条记录, {len(df.columns)} 列")
        return results

    # ── 查询辅助 ──

    def list_files(self, file_type: str) -> list[str]:
        """列出指定类型的所有文件名（已排序）"""
        if file_type not in FILE_TYPES:
            raise ValueError(f"未知文件类型: {file_type}")
        pattern = FILE_TYPES[file_type]["pattern"]
        return sorted(f for f in self._dir_listing if pattern.match(f))

    def summary(self) -> str:
        """打印目录中各类型文件数量"""
        lines = [f"目录: {self.dirpath}"]
        for ftype, spec in FILE_TYPES.items():
            count = sum(1 for f in self._dir_listing if spec["pattern"].match(f))
            lines.append(f"  {ftype:8s}: {count} 个文件")
        return "\n".join(lines)


# ──────────────────────────────────────────────
# 5. 使用示例
# ──────────────────────────────────────────────

if __name__ == "__main__":
    dirpath = "/home/xizhang/MRCDFT-master/Labs/40Ca/40Ca_10_2026-07-14-09-44-24"

    mgr = ElemFileManager(dirpath, do_swap_copy=True, skip_existing=True)

    # 查看概况
    print(mgr.summary())

    # 处理全部 → 两个 DataFrame
    results = mgr.process_all()
    df_kern = results["kern"]
    df_r2   = results["R2_2b"]

    # 也可以单独处理一种
    # df_kern = mgr.process("kern")

    # 保存到 CSV
    # df_kern.to_csv("kern_results.csv", index=False)
    # df_r2.to_csv("r2_results.csv", index=False)

    # 按 beta2 筛选
    # df_kern[(df_kern["beta2_l"] == 0.0) & (df_kern["beta2_r"] == 0.1)]

    print("\n--- kern DataFrame 预览 ---")
    print(df_kern.head())
    print(f"\n列名: {list(df_kern.columns)}")
