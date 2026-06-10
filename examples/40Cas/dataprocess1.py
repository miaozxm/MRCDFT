import sys
import os
# 获取脚本绝对路径（关键！）
current_script_path = os.path.abspath(__file__)
# 直接使用正确的 functionboxes 路径
module_path = "/home/xizhang/functionboxes"

# 验证路径是否存在
if not os.path.exists(module_path):
    print(f"Warning: functionboxes path not found: {module_path}")
    # 尝试其他可能的路径
    alternative_paths = [
        os.path.expanduser("~/functionboxes"),
        os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(current_script_path))), "functionboxes")
    ]
    for alt_path in alternative_paths:
        if os.path.exists(alt_path):
            module_path = alt_path
            print(f"Using alternative path: {module_path}")
            break

# 添加到 sys.path（确保在导入前执行）
sys.path.insert(0, module_path)
from functionboxes.function_file_for_calculation import *
from functionboxes.functionbox import *
# directory = "../MRCDFT/MR_CDFT_f90/datas/40Ca/40Ca11_2026-03-23-11-02-53/output/"
import os, re, shutil
os.chdir('./output')
# 批量处理当前目录下所有 .elem 文件（前缀不固定也支持）
for filename in os.listdir('.'):
    if not os.path.isfile(filename):
        continue
    match = re.search(r'(.*?)([+-])(\d+)([+-])(\d+)_([+-])(\d+)([+-])(\d+)\.elem$', filename)
    if match:
        prefix, symbol1, beta2_1, symbol2, beta3_1, symbol3, beta2_2, symbol4, beta3_2 = match.groups()
        if f"{symbol1}{beta2_1}" != f"{symbol3}{beta2_2}":
            newfilename = prefix + symbol3 + beta2_2 + symbol2 + beta3_1 + '_' + symbol1 + beta2_1 + symbol4 + beta3_2 + '.elem'
            shutil.copy(filename, newfilename)          # 复制原文件内容到交换后的新文件名
            # print(f'交换复制完成: {filename} → {newfilename}')
# kern.1D_eMax08.05.12-010+001_-010+001.elem
# Proj_40Ca_kern.1D_eMax08.05.12-010+001_-010+001.elem
df = main_process_for_ker_data(".")
df.to_csv('output.csv', index=False)
axial_qr=[] ##必选的组态
axial_ql=[]
datar = df
features = ['J','iq2','iq1','bet_2','bet_1','gam_2','gam_1','r2_2','r2_1','N','H']
datar= datar[features]
datar.columns=['JJ','q2','q1','bet_2','bet_1','gam_2','gam_1','r2_2','r2_1','nn','hh']
datar['hh']=datar['hh']*datar['nn']
datar0 = datar[datar['JJ']==0]

# datal= pd.read_csv('Se76GCM.out',sep=r'\s+')
datal = df
features = ['J','iq2','iq1','bet_2','bet_1','gam_2','gam_1','r2_2','r2_1','N','H']
datal= datal[features]
datal.columns=['JJ','q2','q1','bet_2','bet_1','gam_2','gam_1','r2_2','r2_1','nn','hh']
datal['hh']=datal['hh']*datal['nn']
datal0 = datal[datal['JJ']==0]

lab_duple_r,overlap_val_r,labels_r,diagr,dataover_r=create_state2(datar0)
lab_duple_l,overlap_val_l,labels_l,diagl,dataover_l=create_state2(datal0)

split_length=4
Lcutoff=np.arange(0,1.001,0.001)
selr=[]
sell=[]
for i in Lcutoff:
    columns=['bet_2','gam_2','r2_2'] ###形变参数
    selected_r,selected_l,diagdfr, diagdfl=select_Lc4(i,dataover_r,dataover_l,lab_duple_r,
                                     overlap_val_r,lab_duple_l,overlap_val_l,axial_qr,axial_ql,split_length,columns)
    selr.append(selected_r)
    sell.append(selected_l)
# print(diagdfl)
# "            ${mq}             ! number of mesh-point in q-space" "betgam.dat"
selected_l['combined_1'] = selected_l['bet_2'].round(2).astype("str") + "     " + selected_l['gam_2'].round(2).astype("str")
templist = selected_l['combined_1'].tolist()
# 加个元素在templist的最前面

templist.insert(0, f"            {len(templist)}             ! number of mesh-point in q-space")

string = f'\n'.join(templist)
print(string)
# 保存成文件
print(os.getcwd())
# os.mkdir('../exec')
with open("./betgam.dat", "w") as f:
    f.write(string)
    