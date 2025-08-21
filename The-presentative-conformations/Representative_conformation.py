# 轨迹文件和拓扑文件
xtc_file = "step5_1.xtc"
tpr_file = "step5_1.gro"

# 提取帧的时间区间(ps)
period = [
    [2000, 5000,],              # 时间区间，如 2000 ps 到 5000 ps 
    [60000, 120000, 1000],      # 时间区间，如 60000 ps 到 12000 ps，步长为 1000 ps
    80000,                      # 时间点，如 8000 ps
    56000,
    [13560, 35700],
    [47000, 98000],
]

# 结构联配/对齐的时刻（通常选极值点）
reference = 50000
# 最小二乘拟合的参考组 least squares fit group
lsq_fit_group = 4
# 几何中心居中参考组 centering group
centering_group = 1
# 输出轨迹的参考组 output group
extract_group = 1

# Group  0 (      System) has N elements
# Group  1 (     Protein) has N elements
# Group  2 (   Protein-H) has N elements
# Group  3 (     C-alpha) has N elements
# Group  4 (    Backbone) has N elements
# Group  5 (   MainChain) has N elements
# Group  6 (MainChain+Cb) has N elements
# Group  7 ( MainChain+H) has N elements
# Group  8 (   SideChain) has N elements
# Group  9 ( SideChain-H) has N elements
# Group 10 ( Prot-Masses) has N elements
# Group 11 ( non-Protein) has N elements
# Group 12 (       Other) has N elements
# Group 13 (         GTP) has N elements
# Group 14 (          MG) has N elements
# Group 15 (         POT) has N elements
# Group 16 (         CLA) has N elements
# Group 17 (        TIP3) has N elements
# Group 18 (         Ion) has N elements

######################################################
#################### 开发人员分界线 ####################
######################################################

print("\n" + "*"*10 + "\n 开始执行\n" + "*"*10 + "\n")

import shutil
from pathlib import Path
import subprocess
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.rms import rmsd
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")

folder = Path('custom_conformations')
folder.mkdir(exist_ok=True)

print("#####\n# 1 # 提取自定义时段构象\n#####\n")

for frac in period:
    # 执行命令（严格串行）
    try:
        if isinstance(frac, list) and len(frac) == 2:
            outfilename = f"./custom_conformations/_extracted_{frac[0]}_{frac[1]}.xtc"
            # 使用 subprocess.run 替代 Popen 简化同步逻辑
            result = subprocess.run(
                [
                    "gmx", "trjconv",
                    "-s", tpr_file,
                    "-f", xtc_file,
                    "-b", str(frac[0]),
                    "-e", str(frac[1]),
                    "-o", outfilename,
                ],
                input=f"{extract_group}\n".encode(),  # 一次性发送输入
                stdout=subprocess.PIPE,     # 可选：捕获输出
                stderr=subprocess.PIPE,     # 可选：捕获错误
                check=True,               # <=== 加上这个！
            )
            print(result.stdout.decode())
            print(f'[gmx_trjconv] 时段 {frac} {outfilename} 已写入\n'+'-'*10)

        elif isinstance(frac, list) and len(frac) == 3:
            outfilename = f"./custom_conformations/_extracted_{frac[0]}_{frac[1]}_{frac[2]}.xtc"
            # 使用 subprocess.run 替代 Popen 简化同步逻辑
            result = subprocess.run(
                [
                    "gmx", "trjconv",
                    "-s", tpr_file,
                    "-f", xtc_file,
                    "-b", str(frac[0]),
                    "-e", str(frac[1]),
                    "-dt", str(frac[2]),
                    "-o", outfilename,
                ],
                input=f"{extract_group}\n".encode(),  # 一次性发送输入
                stdout=subprocess.PIPE,     # 可选：捕获输出
                stderr=subprocess.PIPE,     # 可选：捕获错误
                check=True,               # <=== 加上这个！
            )
            print(result.stdout.decode())
            print(f'[gmx_trjconv] 时段 {frac}: {outfilename} 已写入\n'+'-'*10)

        elif isinstance(frac, int):
            outfilename = f"./custom_conformations/_extracted_frame_{frac}.xtc"
            # 使用 subprocess.run 替代 Popen 简化同步逻辑
            result = subprocess.run(
                [
                    "gmx", "trjconv",
                    "-s", tpr_file,
                    "-f", xtc_file,
                    "-dump", str(frac),
                    "-o", outfilename,
                ],
                input=f"{extract_group}\n".encode(),  # 一次性发送输入
                stdout=subprocess.PIPE,     # 可选：捕获输出
                stderr=subprocess.PIPE,     # 可选：捕获错误
                check=True,               # <=== 加上这个！
            )
            print(result.stdout.decode())
            print(f'[gmx_trjconv] 时段 {frac} {outfilename} 已写入\n'+'-'*10)
        
        else:
            print(f'[gmx_trjconv] 不认识的时间区段, 跳过\n'+'-'*10)

    except subprocess.CalledProcessError as e:
        print(f"命令执行失败，返回码 {e.returncode}")
        print("错误输出:", e.stderr.decode())
        continue  # 或终止循环（raise）
    except subprocess.TimeoutExpired:
        print("命令执行超时！")
        continue
    except Exception as e:
        print(f"未知异常: {e}")
        continue

# 合并  
pdb_files = list(folder.glob("_extracted*.xtc"))

try:
# 使用 subprocess.run 替代 Popen 简化同步逻辑
    result = subprocess.run(
        [
            "gmx", "trjcat", "-cat",
            "-f", *map(str, pdb_files),
            "-o", '1_1_extracted_combined.xtc',
        ],
        stdout=subprocess.PIPE,     # 可选：捕获输出
        stderr=subprocess.PIPE,     # 可选：捕获错误
        check=True,               # <=== 加上这个！
    )
    print(result.stdout.decode())

    outfilename = Path("./custom_conformations/1_1_extracted_combined.xtc")
    # 如果文件已存在，先删掉
    if outfilename.exists():
        outfilename.unlink()
    shutil.move("1_1_extracted_combined.xtc", str(outfilename))
    print(f'[gmx_trjcat] 已合并所有构象到 {outfilename}\n'+'-'*10)

except subprocess.CalledProcessError as e:
    print("STDOUT:\n", e.stdout)
    print("STDERR:\n", e.stderr)   # ← GROMACS 的真实报错在这里
    raise

for fi in pdb_files:
    fi.unlink()

print(f'[python Path.unlink()] 已删除构象提取临时文件 {pdb_files}\n'+'-'*10)

try:
    outfilename = "./custom_conformations/1_2_reference_frac.pdb"
    # 使用 subprocess.run 替代 Popen 简化同步逻辑
    result = subprocess.run(
        [
            "gmx", "trjconv",
            "-s", tpr_file,
            "-f", xtc_file,
            "-dump", str(reference),
            "-o", outfilename,
        ],
        input=f"{extract_group}\n".encode(),  # 一次性发送输入
        stdout=subprocess.PIPE,     # 可选：捕获输出
        stderr=subprocess.PIPE,     # 可选：捕获错误
        check=True,               # <=== 加上这个！
    )
    print(result.stdout.decode())
    print(f'[gmx_trjconv] 参考匹配构象 {reference}: {outfilename} 已写入\n'+'-'*10)

except subprocess.CalledProcessError as e:
    print("STDOUT:\n", e.stdout)
    print("STDERR:\n", e.stderr)   # ← GROMACS 的真实报错在这里
    raise

print()

############################################################3
print("#####\n# 2 # 帧与参考结构对齐\n#####\n")

# gmx trjconv -f extracted.xtc -s ref.pdb -o aligned.xtc -fit rot+trans -pbc nojump -center

try:
    outfilename = f"./custom_conformations/2_1_aligned.xtc"
    # 使用 subprocess.run 替代 Popen 简化同步逻辑
    result = subprocess.run(
        [
            "gmx", "trjconv",
            "-s", "./custom_conformations/1_2_reference_frac.pdb",
            "-f", "./custom_conformations/1_1_extracted_combined.xtc",
            "-o", outfilename,
            "-fit", "rot+trans", "-center",
        ],
        input=f"{lsq_fit_group}\n{centering_group}\n{extract_group}\n".encode(),  # 一次性发送输入
        stdout=subprocess.PIPE,     # 可选：捕获输出
        stderr=subprocess.PIPE,     # 可选：捕获错误
        check=True,               # <=== 加上这个！
    )
    print(result.stdout.decode())
    print(f'[gmx_trjconv] 各帧与参考结构对齐, {outfilename} 已写入\n'+'-'*10)
        

except subprocess.CalledProcessError as e:
    print(f"命令执行失败，返回码 {e.returncode}")
    print("错误输出:", e.stderr.decode())
except subprocess.TimeoutExpired:
    print("命令执行超时！")
except Exception as e:
    print(f"未知异常: {e}")

print()

######################################################3

print("#####\n# 3 # 计算平均虚拟构象和最近构象\n#####\n")

def average_positions(u, selection="backbone", start=None, stop=None, step=None):
    """
    计算给定原子选择在指定帧窗口上的平均坐标（Å）
    """
    ag = u.select_atoms(selection)
    sum_xyz = np.zeros((ag.n_atoms, 3), dtype=np.float64)
    n = 0
    for ts in u.trajectory[start:stop:step]:
        sum_xyz += ag.positions  # Å
        n += 1
    if n == 0:
        raise ValueError("选定的帧窗口为空，请检查 start/stop/step。")
    return sum_xyz / n, selection

def closest_frame_to_average(u, mean_xyz, selection, start=None, stop=None, step=None, return_frame=False):
    """
    计算每一帧相对平均构象的RMSD（Å），返回最小RMSD对应的帧索引与时间。
    注意：假设 u 对齐过（aligned.xtc），故不再做逐帧拟合。
    """
    ag = u.select_atoms(selection)
    if ag.n_atoms != mean_xyz.shape[0]:
        raise ValueError("mean_xyz 原子数与选择集不匹配。")

    best_idx, best_time, best_rmsd = None, None, np.inf
    rmsd_list = []

    for ts in u.trajectory[start:stop:step]:
        curr = ag.positions  # Å
        val = rmsd(curr, mean_xyz, center=False)  # 已对齐，无需居中/拟合
        rmsd_list.append(val)
        if val < best_rmsd:
            best_rmsd = val
            best_idx = ts.frame
            best_time = ts.time  # ps

    if return_frame:
        return best_idx, best_time, best_rmsd, np.asarray(rmsd_list)
    return best_idx, best_time, best_rmsd

def write_frame_as_pdb(u, frame_index, outfile, selection="all"):
    """
    将指定帧写为PDB（默认写全体系；可用 selection 控制输出原子集）
    """
    u.trajectory[frame_index]
    ag = u.atoms if selection == "all" else u.select_atoms(selection)
    with mda.Writer(outfile, multiframe=False) as W:
        W.write(ag)

# ====== 用法示例 ======
# aligned.xtc 已经对齐；reference_frac.pdb 仅提供拓扑（原子顺序/名字）
u = mda.Universe("./custom_conformations/1_2_reference_frac.pdb", "./custom_conformations/2_1_aligned.xtc")

sel = "protein"      # 或 "name CA" / "protein"
start, stop, step = None, None, None   # 可按需要设置帧窗口与步长

# 🔑 常用选择关键字
# 整体
# "all" → 所有原子
# "protein" → 蛋白质（包括所有原子）
# "nucleic" → 核酸（DNA/RNA）
# "backbone" → 蛋白主链 (N, CA, C, O)
# "name CA" → 仅 Cα 原子
# "resname LIG" → 残基名为 LIG 的小分子/配体
# "segid A" → 段 ID 为 A 的分子

# 逻辑组合
# "protein and name CA" → 蛋白里的 Cα 原子
# "backbone or resname LIG" → 蛋白主链 + LIG 配体
# "protein and not name H*" → 蛋白但不含氢

# 按编号
# "resid 10" → 第 10 号残基
# "resid 10:20" → 残基 10–20
# "bynum 1:1000" → 原子编号 1–1000

# 几何条件
# "around 5 protein" → 蛋白 5 Å 范围内的原子
# "point 10 20 30 5" → 距离点 (10,20,30) 5 Å 内的原子

# 化学类别
# "hydrogen" → 所有氢原子
# "heavy" → 非氢原子
# "polar" / "apolar" → 极性/非极性原子
# "charged" → 带电原子

# 1) 平均构象
mean_xyz, sel_used = average_positions(u, selection=sel, start=start, stop=stop, step=step)

# 2) 逐帧 RMSD 并找最接近的帧
best_idx, best_time_ps, best_rmsd_A, rmsd_series = closest_frame_to_average(
    u, mean_xyz, sel_used, start=start, stop=stop, step=step, return_frame=True,
)

print(f"Closest frame: index={best_idx}, time={best_time_ps:.3f} ps, RMSD={best_rmsd_A:.3f} Å")

# 可选：也把平均构象写出来（只含所选原子集）
sel_ag = u.select_atoms(sel_used)
sel_ag.positions = mean_xyz
with mda.Writer("./custom_conformations/3_1_avg_selected.pdb", multiframe=False) as W:
    W.write(sel_ag)

# 3) 导出最接近帧（整体系或同一选择集）
write_frame_as_pdb(u, best_idx, "./custom_conformations/3_2_closest_to_avg.pdb", selection=sel_used)   # 或 selection=sel_used

print("\n[Python MDAnalysis] 最近构象 ./custom_conformations/3_2_closest_to_avg.pdb 已写入\n")
print('-'*10)
print("\nPowered by Sandy, ChatGPT 5, ChatGPT 5 Thinking\n")
print("Aug 2025\n")
print('-'*10)
print("\n" + "*"*10 + "\n 运行结束 :)\n" + "*"*10 + "\n")

