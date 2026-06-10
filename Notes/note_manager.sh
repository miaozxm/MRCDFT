#!/bin/bash
# ============================================
# MRCDFT 实验笔记管理工具
# ============================================

NOTES_DIR="/home/xizhang/MRCDFT/Notes"
INDEX_FILE="${NOTES_DIR}/index.csv"
TEMPLATES_DIR="${NOTES_DIR}/templates"

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ============================================
# 函数：生成下一个实验ID
# ============================================
generate_exp_id() {
    if [ ! -f "$INDEX_FILE" ]; then
        echo "EXP001"
        return
    fi

    # 从 index.csv 中提取最大的 ID 号
    local max_num=0
    while IFS=',' read -r exp_id rest; do
        if [[ "$exp_id" =~ ^EXP([0-9]+)$ ]]; then
            local num="${BASH_REMATCH[1]}"
            if [ "$num" -gt "$max_num" ]; then
                max_num=$num
            fi
        fi
    done < "$INDEX_FILE"

    local next_num=$((max_num + 1))
    printf "EXP%03d" $next_num
}

# ============================================
# 函数：创建新笔记
# ============================================
create_note() {
    local title="$1"
    local nucleus="$2"
    local nf="$3"
    local template="${4:-default}"
    local group_id="${5:-${nucleus}_study}"

    if [ -z "$title" ] || [ -z "$nucleus" ] || [ -z "$nf" ]; then
        echo -e "${RED}错误: 缺少必要参数${NC}"
        echo "用法: $0 create <title> <nucleus> <nf> [template] [group_id]"
        echo "示例: $0 create 'Ca40形变研究' Ca 10 default Ca40_deformation"
        exit 1
    fi

    local exp_id=$(generate_exp_id)
    local note_file="${NOTES_DIR}/${exp_id}.md"
    local create_time=$(date +"%Y-%m-%d %H:%M:%S")
    local template_file="${TEMPLATES_DIR}/${template}.md"

    # 检查模板是否存在
    if [ ! -f "$template_file" ]; then
        echo -e "${YELLOW}警告: 模板 '${template}' 不存在，使用默认模板${NC}"
        template_file="${TEMPLATES_DIR}/default.md"
    fi

    # 从模板生成笔记
    sed -e "s/{EXP_ID}/$exp_id/g" \
        -e "s/{TITLE}/$title/g" \
        -e "s/{NUCLEUS}/$nucleus/g" \
        -e "s/{NF}/$nf/g" \
        -e "s/{CREATE_TIME}/$create_time/g" \
        -e "s/{GROUP_ID}/$group_id/g" \
        -e "s/{COMPARISON_TYPE}//g" \
        "$template_file" > "$note_file"

    # 更新索引
    echo "${exp_id},${title},${nucleus},${nf},${group_id},${create_time},${create_time},created" >> "$INDEX_FILE"

    echo -e "${GREEN}成功创建实验笔记:${NC}"
    echo -e "  ID: ${BLUE}${exp_id}${NC}"
    echo -e "  文件: ${note_file}"
    echo -e "  请使用此 ID 在 main.sh 中配置: ${YELLOW}NOTE_ID=\"${exp_id}\"${NC}"
}

# ============================================
# 函数：列出所有笔记
# ============================================
list_notes() {
    if [ ! -f "$INDEX_FILE" ]; then
        echo -e "${YELLOW}暂无实验笔记${NC}"
        return
    fi

    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}       实验笔记列表${NC}"
    echo -e "${BLUE}========================================${NC}"
    printf "%-8s %-25s %-8s %-6s %-20s %-10s\n" "ID" "标题" "核素" "Nf" "创建时间" "状态"
    echo -e "${BLUE}----------------------------------------${NC}"

    tail -n +2 "$INDEX_FILE" | while IFS=',' read -r exp_id title nucleus nf group_id created updated status; do
        printf "%-8s %-25s %-8s %-6s %-20s %-10s\n" "$exp_id" "$title" "$nucleus" "$nf" "$created" "$status"
    done

    echo -e "${BLUE}========================================${NC}"
    echo -e "总计: $(tail -n +2 "$INDEX_FILE" | wc -l) 个实验笔记"
}

# ============================================
# 函数：查看指定笔记详情
# ============================================
show_note() {
    local exp_id="$1"

    if [ -z "$exp_id" ]; then
        echo -e "${RED}错误: 请指定实验ID${NC}"
        echo "用法: $0 show <EXP_ID>"
        exit 1
    fi

    local note_file="${NOTES_DIR}/${exp_id}.md"

    if [ ! -f "$note_file" ]; then
        echo -e "${RED}错误: 找不到实验笔记 ${exp_id}${NC}"
        exit 1
    fi

    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}  实验笔记: ${exp_id}${NC}"
    echo -e "${BLUE}========================================${NC}"
    cat "$note_file"
    echo -e "${BLUE}========================================${NC}"
}

# ============================================
# 函数：更新笔记中的实验记录
# ============================================
update_note_record() {
    local exp_id="$1"
    local run_id="$2"
    local run_path="$3"
    local run_time="$4"
    local status="${5:-completed}"

    if [ -z "$exp_id" ] || [ -z "$run_id" ] || [ -z "$run_path" ]; then
        echo -e "${RED}错误: 缺少必要参数${NC}"
        echo "用法: $0 update <EXP_ID> <run_id> <run_path> [run_time] [status]"
        exit 1
    fi

    local note_file="${NOTES_DIR}/${exp_id}.md"

    if [ ! -f "$note_file" ]; then
        echo -e "${RED}错误: 找不到实验笔记 ${exp_id}${NC}"
        exit 1
    fi

    # 在表格中添加新行
    local new_line="| ${run_id} | ${run_path} | ${run_time} | ${status} |"

    # 策略：先删除所有空表格行，然后在表格分隔行后追加新行
    awk -v new_line="$new_line" '
        BEGIN { added = 0 }
        # 跳过空表格行（只包含 | 和空格的行）
        /^\|/ && /\|[[:space:]]*$/ {
            test_line = $0
            gsub(/[|[:space:]]/, "", test_line)
            if (test_line == "") next
        }
        # 打印所有非空行
        { print }
        # 在表格分隔行后添加新行（只添加一次）
        /^ *\| *-.*-.*-.*-.*\|$/ && !added {
            print new_line
            added = 1
        }
    ' "$note_file" > "${note_file}.tmp"
    mv "${note_file}.tmp" "$note_file"

    # 更新索引中的最后更新时间
    local current_time=$(date +"%Y-%m-%d %H:%M:%S")
    if [ -f "$INDEX_FILE" ]; then
        # 使用 awk 替代 sed 以避免路径中的 / 冲突
        awk -F',' -v exp_id="$exp_id" -v time="$current_time" -v status="$status" \
            'BEGIN{OFS=","} $1==exp_id {$7=time; $8=status} {print}' "$INDEX_FILE" > "${INDEX_FILE}.tmp"
        mv "${INDEX_FILE}.tmp" "$INDEX_FILE"
    fi

    echo -e "${GREEN}已更新笔记 ${exp_id}:${NC}"
    echo -e "  添加实验记录: ${run_id}"
}

# ============================================
# 函数：查找同组实验
# ============================================
find_group() {
    local group_id="$1"

    if [ -z "$group_id" ]; then
        echo -e "${RED}错误: 请指定实验组ID${NC}"
        echo "用法: $0 find_group <group_id>"
        exit 1
    fi

    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}  实验组: ${group_id}${NC}"
    echo -e "${BLUE}========================================${NC}"

    if [ -f "$INDEX_FILE" ]; then
        grep ",${group_id}," "$INDEX_FILE" | while IFS=',' read -r exp_id title nucleus nf gid created updated status; do
            echo -e "  ${GREEN}${exp_id}${NC}: ${title} (${nucleus}, Nf=${nf})"
        done
    fi

    echo -e "${BLUE}========================================${NC}"
}

# ============================================
# 主程序
# ============================================
case "$1" in
    create)
        create_note "$2" "$3" "$4" "$5" "$6"
        ;;
    list)
        list_notes
        ;;
    show)
        show_note "$2"
        ;;
    update)
        update_note_record "$2" "$3" "$4" "$5" "$6"
        ;;
    find_group)
        find_group "$2"
        ;;
    *)
        echo -e "${BLUE}========================================${NC}"
        echo -e "${BLUE}  MRCDFT 实验笔记管理工具${NC}"
        echo -e "${BLUE}========================================${NC}"
        echo ""
        echo "用法: $0 <command> [arguments]"
        echo ""
        echo "命令:"
        echo "  create <title> <nucleus> <nf> [template] [group_id]"
        echo "      创建新实验笔记"
        echo "      模板选项: default, comparison"
        echo ""
        echo "  list"
        echo "      列出所有实验笔记"
        echo ""
        echo "  show <EXP_ID>"
        echo "      查看指定实验笔记详情"
        echo ""
        echo "  update <EXP_ID> <run_id> <run_path> [time] [status]"
        echo "      更新实验笔记的记录表格"
        echo ""
        echo "  find_group <group_id>"
        echo "      查找同组的所有实验"
        echo ""
        echo "示例:"
        echo "  $0 create 'Ca40形变研究' Ca 10"
        echo "  $0 create 'Ca40对比实验' Ca 10 comparison Ca40_comparison"
        echo "  $0 list"
        echo "  $0 show EXP001"
        echo "  $0 update EXP001 40Ca_10_xxx /path/to/exp '2026-06-08 17:00' completed"
        echo ""
        ;;
esac
