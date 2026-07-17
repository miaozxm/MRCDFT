#!/bin/bash
# 激活 conda 环境后启动 Qoder IDE
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mrcdft-env
# 启动 IDE（根据实际情况修改命令）
qoder-cn .
