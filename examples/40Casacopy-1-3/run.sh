#!/bin/bash                 		
#PBS -N MRCDFT_Ca_8
#PBS -l select=1:ncpus=11:host=cn2 
#PBS -j oe
                
export OMP_NUM_THREADS=5
export MKL_NUM_THREADS=5
export MKL_DYNAMIC=FALSE
# export MKL_THREADING_LAYER=GNU

# 优先使用PBS环境变量，否则使用脚本真实路径
if [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR" || exit 1
else
    SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
    cd "$SCRIPT_DIR" || exit 1
fi
pwd

# 创建日志目录
LOG_DIR="$(pwd)/logs"
mkdir -p "$LOG_DIR"

# 获取当前时间戳用于日志文件名
current_time=$(date +"%Y%m%d_%H%M%S")

# 定义日志文件路径
LOG_FILE="$LOG_DIR/run_${current_time}.log"
ERROR_FILE="$LOG_DIR/error_${current_time}.log"

# 开始重定向输出
# exec > >(tee -a "$LOG_FILE") 2>&1
# 或者分开重定向
exec 1> >(tee -a "$LOG_FILE")
exec 2> >(tee -a "$ERROR_FILE")

echo "日志文件: $LOG_FILE"
echo "开始时间: $(date)"

start_time=$(date +%s)
echo -e "\033[32m run ...\033[0m"

# 使用MRCDFT的绝对路径
MRCDFT_BIN="/home/xizhang/MRCDFT/bin/MRCDFT"
mpirun -np 15 "$MRCDFT_BIN" -p para.dat -d b23.dat

echo calculation is finished !
end_time=$(date +%s)

execution_time=$((end_time - start_time))
execution_time_minutes=$((execution_time / 60))
execution_time_seconds=$((execution_time % 60))
echo "Time cost : ${execution_time_minutes}min${execution_time_seconds}s"

echo Done!



