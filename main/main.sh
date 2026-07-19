#!/bin/bash                 		
#PBS -N MRCDFT
#PBS -l select=1:ncpus=63
#PBS -o /dev/null
#PBS -e /dev/null

FDIR="/home/xizhang/MRCDFT-master"
# ===== 基础路径配置（唯一绝对路径）=====
MAIN_DIR="$FDIR/main"
LABS_DIR="$FDIR/Labs"
MRCDFT_BIN="$FDIR/bin/MRCDFT"

cd "$MAIN_DIR" || { echo "ERROR: Cannot cd to $MAIN_DIR"; exit 1; }

export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export MKL_DYNAMIC=FALSE
# export MKL_THREADING_LAYER=GNU

# ===== 实验参数配置 =====
ELE="Ca"
A=48

# ===== 实验笔记ID配置（可选）=====
# 如果设置了 NOTE_ID，运行结束后会自动更新对应笔记
# 留空则不更新笔记系统
# EXP_PURPOSE="${1:-}"
EXP_PURPOSE="r2-2body"
NOTE_ID=""  # 例如: "EXP001"

# 从 para.dat 自动提取 Nf
Nf=$(grep -E '^\s*n0f' para.dat | head -1 | awk '{print $3}')
if [ -z "$Nf" ]; then
    echo "ERROR: Cannot extract Nf from para.dat"
    exit 1
fi

# ===== 同步更新 para.dat 中的核素信息 =====
sed -i "s/^[[:space:]]*[A-Za-z][a-z]*[[:space:]]*[0-9]\+[[:space:]]*!.*nucleus/${ELE} ${A}                               ! nucleus/" para.dat

echo "已更新 para.dat: ${ELE} ${A}"

# ===== 创建实验空间 =====
current_time=$(date +"%Y-%m-%d-%H-%M-%S")
log_time=$current_time
NUC="${A}${ELE}"
EXP_ID="${NUC}_${Nf}_${current_time}"
if [ -n "$EXP_PURPOSE" ]; then
    EXP_PATH="${LABS_DIR}/${NUC}/${EXP_PURPOSE}/${EXP_ID}"
else
    EXP_PATH="${LABS_DIR}/${NUC}/${EXP_ID}"
fi

# 创建目录结构
mkdir -p "${EXP_PATH}/output"
mkdir -p "${EXP_PATH}/logs"

# 复制配置文件（此时 para.dat 已经是更新后的版本）
cp para.dat "${EXP_PATH}/"
cp b23.dat "${EXP_PATH}/"
[ -f para2.dat ] && cp para2.dat "${EXP_PATH}/"

# ===== 实验笔记管理 =====
NOTE_FILE="experiment_note.txt"
if [ -f "$NOTE_FILE" ] && [ -s "$NOTE_FILE" ]; then
    # 如果笔记文件存在且非空，直接复制
    cp "$NOTE_FILE" "${EXP_PATH}/"
    echo "已复制实验笔记: $NOTE_FILE"
else
    # 否则生成默认笔记
    cat > "${EXP_PATH}/${NOTE_FILE}" <<EOF
# 实验记录
# ============================================
核素: ${NUC}
谐振子壳层数 (Nf): ${Nf}
实验ID: ${EXP_ID}
运行时间: $(date +"%Y-%m-%d %H:%M:%S")
============================================

## 实验目的
老src测试352.99来源

## 相关实验
- 参考实验：（如有，请填写实验ID）
- 对比实验：（如有，请填写实验ID）

## 备注
（其他需要记录的信息）
EOF
    echo "已生成默认实验笔记: ${NOTE_FILE}"
fi

# 保存实验组ID（用于关联同系列实验）
GROUP_ID="${NUC}_study"
echo "$GROUP_ID" > "${EXP_PATH}/.group_id"

# 切换到实验目录
cd "${EXP_PATH}" || { echo "ERROR: Cannot cd to $EXP_PATH"; exit 1; }
echo "===== 实验目录: ${EXP_PATH} ====="

# ===== 日志设置 =====
LOG_DIR="./logs"
# echo $log_time
LOG_FILE="$LOG_DIR/run_${log_time}.log"
ERROR_FILE="$LOG_DIR/error_${log_time}.log"

exec 1> >(tee -a "$LOG_FILE")
exec 2> >(tee -a "$ERROR_FILE")

# 将当前脚本复制到logs目录，并修改路径参数使其可直接在目标位置重复运行
# 1. MAIN_DIR 改为 ..（实验目录，参数文件所在位置）
# 2. current_time 固定为原始运行时间（防止目录重定向到新时间戳）
sed -e "s|^MAIN_DIR=.*|MAIN_DIR=\"..\"|" \
    -e "s|^current_time=.*|current_time=\"${current_time}\"|" \
    "$0" > "$LOG_DIR/run_${log_time}.sh"
chmod +x "$LOG_DIR/run_${log_time}.sh"

echo "日志文件: $LOG_FILE"
echo "开始时间: $(date)"

start_time=$(date +%s)
echo -e "\033[32m run ...\033[0m"

mpirun -np 21 "$MRCDFT_BIN" -p para.dat -d b23.dat
# mpirun -np 28 "$MRCDFT_BIN" -p para2.dat -d b23.dat

echo calculation is finished !
end_time=$(date +%s)

execution_time=$((end_time - start_time))
execution_time_minutes=$((execution_time / 60))
execution_time_seconds=$((execution_time % 60))
echo "Time cost : ${execution_time_minutes}min${execution_time_seconds}s"

echo Done!

# ===== 实验后处理：更新笔记系统 =====
if [ -n "$NOTE_ID" ]; then
    NOTE_MANAGER="/home/xizhang/MRCDFT/Notes/note_manager.sh"
    if [ -x "$NOTE_MANAGER" ]; then
        run_time=$(date + "%Y-%m-%d %H:%M:%S")
        echo "正在更新实验笔记 ${NOTE_ID}..."
        bash "$NOTE_MANAGER" update "$NOTE_ID" "$EXP_ID" "$EXP_PATH" "$run_time" "completed"
    else
        echo "警告: 找不到笔记管理工具，跳过笔记更新"
    fi
fi

# ===== 实验后清理：重置笔记文件 =====
# 如果 main/experiment_note.txt 存在，清空它以触发下次运行时的默认生成
if [ -f "${MAIN_DIR}/${NOTE_FILE}" ]; then
    > "${MAIN_DIR}/${NOTE_FILE}"
    echo "已重置实验笔记模板: ${NOTE_FILE}"
fi



