# MRCDFT 实验笔记系统

## 概述

独立的实验笔记管理系统，用于追踪和管理 MRCDFT 计算实验。

## 目录结构

```
Notes/
├── note_manager.sh      # 笔记管理工具
├── index.csv            # 笔记索引文件
├── templates/           # 笔记模板
│   ├── default.md       # 默认模板
│   └── comparison.md    # 对比实验模板
├── EXP001.md            # 实验笔记（自动生成）
├── EXP002.md
└── README.md            # 本文件
```

## 快速开始

### 1. 创建实验笔记

```bash
# 基本用法
./note_manager.sh create "实验标题" 核素 Nf

# 示例：创建默认笔记
./note_manager.sh create "Ca40形变能量曲线研究" Ca 10

# 使用对比实验模板
./note_manager.sh create "Ca40对比实验" Ca 10 comparison Ca40_comparison
```

### 2. 在 main.sh 中配置笔记 ID

编辑 `main/main.sh`，设置 `NOTE_ID`：

```bash
NOTE_ID="EXP001"  # 填入创建的笔记ID
```

### 3. 运行实验

```bash
cd main
bash main.sh
# 或
qsub main.sh
```

运行结束后，脚本会自动更新笔记中的实验记录表格。

## 命令参考

### create - 创建新笔记

```bash
./note_manager.sh create <title> <nucleus> <nf> [template] [group_id]
```

**参数：**
- `title`: 实验标题
- `nucleus`: 核素符号（如 Ca）
- `nf`: 谐振子壳层数
- `template`: 模板名称（default 或 comparison，默认为 default）
- `group_id`: 实验组ID（可选，默认为 `{nucleus}_study`）

### list - 列出所有笔记

```bash
./note_manager.sh list
```

### show - 查看笔记详情

```bash
./note_manager.sh show <EXP_ID>
# 示例
./note_manager.sh show EXP001
```

### update - 更新笔记记录

```bash
./note_manager.sh update <EXP_ID> <run_id> <run_path> [time] [status]
```

**注意：** 此命令通常由 main.sh 自动调用，无需手动执行。

### find_group - 查找同组实验

```bash
./note_manager.sh find_group <group_id>
# 示例
./note_manager.sh find_group Ca40_deformation
```

## 工作流程

### 场景 1：单次实验

```bash
# 1. 创建笔记
./note_manager.sh create "Ca40基准计算" Ca 10

# 2. 在 main.sh 中设置 NOTE_ID="EXP001"

# 3. 运行实验
bash main.sh

# 4. 笔记自动更新，查看结果
./note_manager.sh show EXP001
```

### 场景 2：对比实验

```bash
# 1. 创建对比实验笔记
./note_manager.sh create "Ca40_Nf对比" Ca 10 comparison Ca40_Nf_comparison

# 2. 第一次运行（Nf=8）
# 在 main.sh 中设置 NOTE_ID="EXP00X", Nf=8
bash main.sh

# 3. 第二次运行（Nf=10）
# 在 main.sh 中设置 NOTE_ID="EXP00X", Nf=10
bash main.sh

# 4. 查看对比记录
./note_manager.sh show EXP00X
```

### 场景 3：查找相关实验

```bash
# 查找所有 Ca40 形变研究的实验
./note_manager.sh find_group Ca40_deformation

# 查看所有笔记
./note_manager.sh list
```

## 笔记模板

### 默认模板 (default.md)

适用于常规实验，包含：
- 基本信息
- 实验目的
- 相关实验
- 实验记录表格
- 结果分析
- 备注

### 对比实验模板 (comparison.md)

适用于对比实验，包含：
- 对比的实验列表
- 实验组 A/B 的记录表格
- 对比结果
- 结论

## 索引文件 (index.csv)

格式：
```csv
EXP_ID,Title,Nucleus,Nf,GroupID,CreatedTime,LastUpdated,Status
```

每次创建或更新笔记时自动维护。

## 最佳实践

1. **为每个研究课题创建独立的笔记**
   - 例如：形变研究、约束计算、投影计算等分别创建笔记

2. **使用有意义的实验组ID**
   - 例如：`Ca40_deformation`、`Ca40_PNP_comparison`

3. **在实验前填写实验目的**
   - 运行前编辑笔记文件，填写实验目的和相关实验

4. **定期查看和整理笔记**
   - 使用 `list` 命令查看所有实验
   - 使用 `find_group` 查找相关实验

5. **备份 Notes 目录**
   - Notes 目录包含所有实验元数据，建议定期备份

## 与旧系统的兼容

本系统与原有的 `experiment_note.txt` 机制兼容：
- 如果设置了 `NOTE_ID`，优先使用笔记系统
- 如果未设置 `NOTE_ID`，回退到原来的 `experiment_note.txt` 机制
