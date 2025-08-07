
# %% 
import pandas as pd
import glob
import matplotlib.pyplot as plt
import YLutilpy
import numpy as np
YLutilpy.default_img_set()
# %% 

# # 读取CSV文件，第一行作为表头
# subID = 1
# groupID = 1
# # 构建文件名前缀
# prefix = f"../Data/V_Result_G{groupID}_Sub{subID}_"
# files = glob.glob(prefix + "*.csv")
# if not files:
#     raise FileNotFoundError(f"No file found with prefix {prefix}")
# filename = files[0]
# df = pd.read_csv(filename)

# # 剔除ID列小于等于0的行
# df = df[df['ID'] > 0]
# 支持多个subID和两种类型（V和A）的文件读取与合并
subIDs = [1]  # 这里填你要分析的被试编号
groupID = 1

all_dfs = []
for subID in subIDs:
    for prefix_type in ['V', 'A']:
        prefix = f"../Data/{prefix_type}_Result_G{groupID}_Sub{subID}_"
        files = glob.glob(prefix + "*.csv")
        if not files:
            print(f"Warning: No file found with prefix {prefix}")
            continue
        filename = files[0]
        df_tmp = pd.read_csv(filename)
        df_tmp = df_tmp[df_tmp['ID'] > 0]
        df_tmp['subID'] = subID
        df_tmp['modality'] = prefix_type
        all_dfs.append(df_tmp)

if not all_dfs:
    raise FileNotFoundError("No data files found for any subject or modality.")

df = pd.concat(all_dfs, ignore_index=True)

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
fig, axes = plt.subplots(1, 3, figsize=(12, 5), dpi=300)

# RT柱状图
pivot_rt.plot(kind='bar', ax=axes[0])
axes[0].set_xlabel('tSOA')
axes[0].set_ylabel('RT')
axes[0].set_title('Mean RT by cueType and tSOA')
# axes[0].legend(title='cueType')

# judge柱状图
pivot_judge.plot(kind='bar', ax=axes[1])
axes[1].set_xlabel('tSOA')
axes[1].set_ylabel('judge')
axes[1].set_title('ACC by cueType and tSOA')
# axes[1].legend(title='cueType')

# RT/judge比值柱状图
pivot_ratio = pivot_rt / pivot_judge
ax=axes[2]
pivot_ratio.plot(kind='bar', ax=ax)
ax.set_xlabel('tSOA')
ax.set_ylabel('RT / ACC')
ax.set_title('RT/ACC Ratio by cueType and tSOA')
ax.legend(title='cueType',loc='upper left', bbox_to_anchor=(1, 1))
axes[0].get_legend().remove()
axes[1].get_legend().remove()

plt.tight_layout()
plt.show()
# %%
# 自定义颜色
color_dict = {
    'AP1': '#FFE699',   
    'AP2': '#FFD966',   
    'PP':  '#42A62A',   
    'AU':  "#EC7025",   
}

# 应用到前面所有绘图
# 1. summary 柱状图
cue_order = ['AP1', 'AP2', 'PP', 'AU']
summary = summary.set_index('cueType').reindex(cue_order).reset_index()
plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
plt.bar(summary['cueType'], summary['RT'], color=[color_dict.get(c, 'gray') for c in summary['cueType']])
plt.xlabel('cueType')
plt.ylabel('RT')
plt.title('Mean RT by cueType')

plt.subplot(1,2,2)
plt.bar(summary['cueType'], summary['judge'], color=[color_dict.get(c, 'gray') for c in summary['cueType']])
plt.xlabel('cueType')
plt.ylabel('judge')
plt.title('Mean judge by cueType')

plt.tight_layout()
plt.show()

# 2. summary_tSOA 折线图
fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)
for cue in cue_order:
    data = summary_tSOA[summary_tSOA['cueType'] == cue]
    axes[0].plot(data['tSOA'], data['RT'], marker='o', label=f'cueType {cue}', color=color_dict.get(cue, 'gray'))
axes[0].set_xlabel('tSOA')
axes[0].set_ylabel('RT')
axes[0].set_title('Mean RT by cueType and tSOA')
axes[0].legend()

for cue in cue_order:
    data = summary_tSOA[summary_tSOA['cueType'] == cue]
    axes[1].plot(data['tSOA'], data['judge'], marker='o', label=f'cueType {cue}', color=color_dict.get(cue, 'gray'))
axes[1].set_xlabel('tSOA')
axes[1].set_ylabel('ACC')
axes[1].set_title('ACC by cueType and tSOA')
axes[1].legend()

for cue in cue_order:
    data = summary_tSOA[summary_tSOA['cueType'] == cue].copy()
    data['RT_judge_ratio'] = data['RT'] / data['judge']
    axes[2].plot(data['tSOA'], data['RT_judge_ratio'], marker='o', label=f'cueType {cue}', color=color_dict.get(cue, 'gray'))
axes[2].set_xlabel('tSOA')
axes[2].set_ylabel('RT / ACC')
axes[2].set_title('RT/ACC Ratio by cueType and tSOA')
axes[2].legend()

plt.tight_layout()
plt.show()

# 3. 分组柱状图
fig, axes = plt.subplots(1, 3, figsize=(12, 5), dpi=300)
pivot_rt = pivot_rt[cue_order]
pivot_judge = pivot_judge[cue_order]
pivot_ratio = pivot_rt / pivot_judge

pivot_rt.plot(kind='bar', ax=axes[0], color=[color_dict.get(c, 'gray') for c in cue_order])
axes[0].set_xlabel('tSOA')
axes[0].set_ylabel('RT')
axes[0].set_title('Mean RT by cueType and tSOA')

pivot_judge.plot(kind='bar', ax=axes[1], color=[color_dict.get(c, 'gray') for c in cue_order])
axes[1].set_xlabel('tSOA')
axes[1].set_ylabel('judge')
axes[1].set_title('ACC by cueType and tSOA')

pivot_ratio.plot(kind='bar', ax=axes[2], color=[color_dict.get(c, 'gray') for c in cue_order])
axes[2].set_xlabel('tSOA')
axes[2].set_ylabel('RT / ACC')
axes[2].set_title('RT/ACC Ratio by cueType and tSOA')
axes[2].legend(title='cueType',loc='upper left', bbox_to_anchor=(1, 1))
axes[0].get_legend().remove()
axes[1].get_legend().remove()

plt.tight_layout()
plt.show()
# %%
# 调整分组柱状图的cueType顺序，使PP在每组最左边
# 只需调整cue_order，把'PP'放在第一个
cue_order_new = ['PP', 'AP1', 'AP2', 'AU']

# 重新排序pivot表
pivot_rt_new = pivot_rt[cue_order_new]
pivot_judge_new = pivot_judge[cue_order_new]
pivot_ratio_new = pivot_rt_new / pivot_judge_new

fig, axes = plt.subplots(1, 3, figsize=(12, 5), dpi=300)
pivot_rt_new.plot(kind='bar', ax=axes[0], color=[color_dict.get(c, 'gray') for c in cue_order_new])
axes[0].set_xlabel('tSOA')
axes[0].set_ylabel('RT')
axes[0].set_title('Mean RT by cueType and tSOA')

pivot_judge_new.plot(kind='bar', ax=axes[1], color=[color_dict.get(c, 'gray') for c in cue_order_new])
axes[1].set_xlabel('tSOA')
axes[1].set_ylabel('judge')
axes[1].set_title('ACC by cueType and tSOA')

pivot_ratio_new.plot(kind='bar', ax=axes[2], color=[color_dict.get(c, 'gray') for c in cue_order_new])
axes[2].set_xlabel('tSOA')
axes[2].set_ylabel('RT / ACC')
axes[2].set_title('RT/ACC Ratio by cueType and tSOA')
axes[2].legend(title='cueType',loc='upper left', bbox_to_anchor=(1, 1))
axes[0].get_legend().remove()
axes[1].get_legend().remove()

plt.tight_layout()
plt.show()

# %% 
# 只保留最后重排cueType顺序并绘图的部分
# 
import matplotlib.pyplot as plt
YLutilpy.default_img_set()

subIDs = [1]  # 这里填你要分析的被试编号
groupID = 1

all_dfs = []
for subID in subIDs:
    for prefix_type in ['V']:
        prefix = f"../Data/{prefix_type}_Result_G{groupID}_Sub{subID}_"
        files = glob.glob(prefix + "*.csv")
        if not files:
            continue
        filename = files[0]
        df_tmp = pd.read_csv(filename)
        df_tmp = df_tmp[df_tmp['ID'] > 0]
        df_tmp['subID'] = subID
        df_tmp['modality'] = prefix_type
        all_dfs.append(df_tmp)

if not all_dfs:
    raise FileNotFoundError("No data files found for any subject or modality.")

df = pd.concat(all_dfs, ignore_index=True)

# 将tSOA分为小、中、大三组
tsoa_bins = pd.qcut(df['tSOA'], 3, labels=['Short Invalid', 'Valid', 'Long Invalid'])
df['tSOA_group'] = tsoa_bins

# 按cueType和tSOA_group分组，计算均值
grouped = df.groupby(['tSOA_group', 'cueType'])[['RT', 'judge']].mean().reset_index()

# 透视表用于绘图
pivot_rt = grouped.pivot(index='tSOA_group', columns='cueType', values='RT')
pivot_judge = grouped.pivot(index='tSOA_group', columns='cueType', values='judge')

# cueType顺序调整，PP最左
cue_order_new = [ 'AP2','PP', 'AP1', 'AU']
color_dict = {
    'AP1': '#FFE699',
    'AP2': '#FFD966',
    'PP':  '#42A62A',
    'AU':  "#EC7025",
}

pivot_rt_new = pivot_rt[cue_order_new]
pivot_judge_new = pivot_judge[cue_order_new]
pivot_ratio_new = pivot_rt_new / pivot_judge_new

fig, axes = plt.subplots(1, 3, figsize=(12, 5), dpi=300)
pivot_rt_new.plot(kind='bar', ax=axes[0], color=[color_dict.get(c, 'gray') for c in cue_order_new])
axes[0].set_xlabel('tSOA')
axes[0].set_ylabel('RT')
axes[0].set_title('Mean RT by cueType and tSOA')

pivot_judge_new.plot(kind='bar', ax=axes[1], color=[color_dict.get(c, 'gray') for c in cue_order_new])
axes[1].set_xlabel('tSOA')
axes[1].set_ylabel('judge')
axes[1].set_title('ACC by cueType and tSOA')

pivot_ratio_new.plot(kind='bar', ax=axes[2], color=[color_dict.get(c, 'gray') for c in cue_order_new])
axes[2].set_xlabel('tSOA')
axes[2].set_ylabel('RT / ACC')
axes[2].set_title('RT/ACC Ratio by cueType and tSOA')
axes[2].legend(title='cueType',loc='upper left', bbox_to_anchor=(1, 1))
axes[0].get_legend().remove()
axes[1].get_legend().remove()

plt.tight_layout()
plt.show()
# %%
# 绘制ID<0的trial中tgAmp的变化曲线
for subID in [1,2,3,4,5]:
    for prefix_type in ['V', 'A']:
        prefix = f"../Data/{prefix_type}_Result_G{groupID}_Sub{subID}_"
        files = glob.glob(prefix + "*.csv")
        if not files:
            print(f"Warning: No file found with prefix {prefix}")
            continue
        filename = files[0]
        df = pd.read_csv(filename)

        df_sub5_invalid = df[(df['ID'] < 0)]
        if df_sub5_invalid.empty:
            print(f"No invalid trials (ID<0) found for subID={subID}.")
        else:
            plt.figure(figsize=(8,4))
            plt.plot(df_sub5_invalid['tgAmp'].values, marker='o')
            plt.xlabel('Trial Index')
            plt.ylabel('tgAmp')
            plt.title(f'tgAmp of subID={subID}, Mod={prefix_type}')
            plt.tight_layout()
            plt.show()
# %%
