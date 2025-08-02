
# %% 
import pandas as pd
import glob
import matplotlib.pyplot as plt
import YLutilpy
import numpy as np
YLutilpy.default_img_set()
# %% 

# 读取CSV文件，第一行作为表头
subID = 1
groupID = 1
# 构建文件名前缀
prefix = f"../Data/V_Result_G{groupID}_Sub{subID}_"
files = glob.glob(prefix + "*.csv")
if not files:
    raise FileNotFoundError(f"No file found with prefix {prefix}")
filename = files[0]
df = pd.read_csv(filename)

# 剔除ID列小于等于0的行
df = df[df['ID'] > 0]


# %% 
# 按cueType分组，统计RT和judge的均值
summary = df.groupby('cueType')[['RT', 'judge']].mean().reset_index()
print(summary)

# 可视化RT
plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
plt.bar(summary['cueType'], summary['RT'], color='skyblue')
plt.xlabel('cueType')
plt.ylabel('RT')
plt.title('Mean RT by cueType')

# 可视化judge
plt.subplot(1,2,2)
plt.bar(summary['cueType'], summary['judge'], color='salmon')
plt.xlabel('cueType')
plt.ylabel('judge')
plt.title('Mean judge by cueType')

plt.tight_layout()
plt.show()

# %%
summary_tSOA = df.groupby(['cueType', 'tSOA'])[['RT', 'judge']].mean().reset_index()

# 按cueType和tSOA分组，统计RT和judge的均值
fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)

# RT by cueType and tSOA
for cue in summary_tSOA['cueType'].unique():
    data = summary_tSOA[summary_tSOA['cueType'] == cue]
    axes[0].plot(data['tSOA'], data['RT'], marker='o', label=f'cueType {cue}')
axes[0].set_xlabel('tSOA')
axes[0].set_ylabel('RT')
axes[0].set_title('Mean RT by cueType and tSOA')
axes[0].legend()

# judge by cueType and tSOA
for cue in summary_tSOA['cueType'].unique():
    data = summary_tSOA[summary_tSOA['cueType'] == cue]
    axes[1].plot(data['tSOA'], data['judge'], marker='o', label=f'cueType {cue}')
axes[1].set_xlabel('tSOA')
axes[1].set_ylabel('ACC')
axes[1].set_title('ACC by cueType and tSOA')
axes[1].legend()

# RT/judge比值 by cueType and tSOA
for cue in summary_tSOA['cueType'].unique():
    data = summary_tSOA[summary_tSOA['cueType'] == cue].copy()
    data['RT_judge_ratio'] = data['RT'] / data['judge']
    axes[2].plot(data['tSOA'], data['RT_judge_ratio'], marker='o', label=f'cueType {cue}')
axes[2].set_xlabel('tSOA')
axes[2].set_ylabel('RT / ACC')
axes[2].set_title('RT/ACC Ratio by cueType and tSOA')
axes[2].legend()

plt.tight_layout()
plt.show()
# %% 
# 将tSOA分为小、中、大三组
tsoa_bins = pd.qcut(df['tSOA'], 3, labels=['Short Invalid', 'Valid', 'Long Invalid'])
df['tSOA_group'] = tsoa_bins

# 按cueType和tSOA_group分组，计算均值
grouped = df.groupby(['tSOA_group', 'cueType'])[['RT', 'judge']].mean().reset_index()

# 透视表用于绘图
pivot_rt = grouped.pivot(index='tSOA_group', columns='cueType', values='RT')
pivot_judge = grouped.pivot(index='tSOA_group', columns='cueType', values='judge')

# 绘图
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

# RT柱状图
pivot_rt.plot(kind='bar', ax=axes[0])
axes[0].set_xlabel('tSOA')
axes[0].set_ylabel('RT')
axes[0].set_title('Mean RT by cueType and tSOA')
axes[0].legend(title='cueType')

# judge柱状图
pivot_judge.plot(kind='bar', ax=axes[1])
axes[1].set_xlabel('tSOA分组')
axes[1].set_ylabel('judge')
axes[1].set_title('ACC by cueType and tSOA')
axes[1].legend(title='cueType')

plt.tight_layout()
plt.show()

# %%
