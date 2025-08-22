
# %% 
import pandas as pd
import glob
import matplotlib.pyplot as plt
import YLutilpy
import numpy as np
from scipy.stats import ttest_rel

YLutilpy.default_img_set()


for subIDs in [[1],[2]]:  # 被试编号
# for subIDs in [[4]]:  # 被试编号
    for prefix_types in [['A'],['V']]:  # 组编号
        groupID = 1
        all_dfs = []

        # 加载数据
        # 初始化存储每个被试的平均值
        rt_means = {}
        acc_means = {}
        ratio_means = {}
        for prefix_type in prefix_types:  # 可扩展到其他模态如'A'
            for subID in subIDs:
                prefix = f"../Data/{prefix_type}_Result_G{groupID}_Sub{subID}_"
                files = glob.glob(prefix + "*.csv")
                if not files:
                    continue
                filename = files[0]
                df_tmp = pd.read_csv(filename)
                df_tmp = df_tmp[df_tmp['ID'] > 0 ]  # 保留有效ID
                df_tmp['subID'] = subID
                df_tmp['modality'] = prefix_type
                # 去除RT 大于或小于2.5个标准差的行
                df_tmp = df_tmp[(df_tmp['RT'] < df_tmp['RT'].mean() + 2.5 * df_tmp['RT'].std()) & 
                                (df_tmp['RT'] > df_tmp['RT'].mean() - 2.5 * df_tmp['RT'].std()) |
                                (df_tmp['RT'].isna())]  # 保留RT为NaN的行（错误反应）
                all_dfs.append(df_tmp)
        if not all_dfs:
            raise FileNotFoundError("No data files found.")

        df = pd.concat(all_dfs, ignore_index=True)

        df['tSOA_group'] = df.groupby('cueType')['tSOA'].rank(method='dense').astype(int)    
        # 将tSOA_group映射成字符标签
        tsoa_group_mapping = {1: 'Short ', 2: 'Long'}
        # 根据cueType和tSOA_group的关系确定validity
        df['validity'] = np.where(
            ((df['cueType'] == 'AUl') & (df['tSOA_group'] == 'Long')), 'valid',
            np.where(
                ((df['cueType'] == 'AUs') & (df['tSOA_group'] == 'Short')), 'invalid',
                np.where(
                    ((df['cueType'] == 'AUs') & (df['tSOA_group'] == 'Short')), 'valid',
                    np.where(
                        ((df['cueType'] == 'AUs') & (df['tSOA_group'] == 'Long')), 'invalid',
                        None
                    )
                )
            )
        )

        # 将数据转换为表格形式
        table_RT = df[df['judge'] == 1].groupby(['subID', 'cueType', 'tSOA_group'])['RT'].mean().unstack(['tSOA_group'])
        table_Acc = df.groupby(['subID', 'cueType', 'tSOA_group'])['judge'].mean().unstack(['tSOA_group'])
        table_ratio = np.divide(table_RT,table_Acc)

# %% 
import pandas as pd
import glob
import matplotlib.pyplot as plt
import YLutilpy
import numpy as np
from scipy.stats import ttest_rel

YLutilpy.default_img_set()

# 支持多个subID和两种类型（V和A）的文件读取与合并[1],[2],[3],[4],[5],
for subIDs in [[2,3,4,5,6,7,8,9,10,11]]:  # 被试编号
# for subIDs in [[4]]:  # 被试编号
    for prefix_types in [['A'],['V']]:  # 组编号
        groupID = 1
        all_dfs = []

        # 加载数据
        # 初始化存储每个被试的平均值
        rt_means = {}
        acc_means = {}
        ratio_means = {}
        for prefix_type in prefix_types:  # 可扩展到其他模态如'A'
            for subID in subIDs:
                prefix = f"../Data_0813/{prefix_type}_Result_G{groupID}_Sub{subID}_"
                files = glob.glob(prefix + "*.csv")
                if not files:
                    continue
                filename = files[0]
                df_tmp = pd.read_csv(filename)
                df_tmp = df_tmp[(df_tmp['ID'] > 0 )* df_tmp['cueType']=='PP']  # 保留有效ID
                df_tmp['subID'] = subID
                df_tmp['modality'] = prefix_type
                # 去除RT 大于或小于2.5个标准差的行
                df_tmp = df_tmp[(df_tmp['RT'] < df_tmp['RT'].mean() + 2.5 * df_tmp['RT'].std()) & 
                                (df_tmp['RT'] > df_tmp['RT'].mean() - 2.5 * df_tmp['RT'].std()) |
                                (df_tmp['RT'].isna())]  # 保留RT为NaN的行（错误反应）
                all_dfs.append(df_tmp)
        if not all_dfs:
            raise FileNotFoundError("No data files found.")

        df = pd.concat(all_dfs, ignore_index=True)

        df['tSOA_group'] = df.groupby('cueType')['tSOA'].rank(method='dense').astype(int)    
        # 将tSOA_group映射成字符标签
        tsoa_group_mapping = {1: 'Short Invalid', 2: 'Valid', 3: 'Long Invalid'}
        df['tSOA_group'] = df['tSOA_group'].map(tsoa_group_mapping)   

        # 将数据转换为表格形式
        table_RT = df[df['judge'] == 1].groupby(['subID', 'tSOA_group'])['RT'].mean().unstack(['tSOA_group'])
        table_Acc = df.groupby(['subID', 'tSOA_group'])['judge'].mean().unstack(['tSOA_group'])
        table_ratio = np.divide(table_RT,table_Acc)

        # 绘柱状图
        fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=300)
        # 调整tSOA_group的顺序为 Short Invalid -> Valid -> Long Invalid
        tsoa_order = ['Short Invalid', 'Valid', 'Long Invalid']
        # 重新排序table column
        table_RT = table_RT.reindex(tsoa_order, axis=1)
        table_Acc = table_Acc.reindex(tsoa_order, axis=1)
        # 对每个表的所有行取平均值
        table_RT_mean = table_RT.mean(axis=0)
        table_Acc_mean = table_Acc.mean(axis=0)
        table_ratio_mean = table_ratio.mean(axis=0)

        # 将平均值转换为DataFrame以便绘图
        table_RT_mean = pd.DataFrame(table_RT_mean).T
        table_Acc_mean = pd.DataFrame(table_Acc_mean).T
        table_ratio_mean = pd.DataFrame(table_ratio_mean).T
        # 绘制RT子图
        table_RT_mean.T.plot(kind='bar', ax=axes[0])
        axes[0].set_title('Mean Reaction Time (Correct Trials Only)')
        axes[0].set_xlabel('tSOA Group')
        axes[0].set_ylabel('RT (s)')
        axes[0].get_legend().remove()

        # 绘制ACC子图
        table_Acc_mean.T.plot(kind='bar', ax=axes[1])
        axes[1].set_title('Accuracy (All Trials)')
        axes[1].set_xlabel('tSOA Group')
        axes[1].set_ylabel('Accuracy')
        axes[1].get_legend().remove()

        # 绘制比值子图
        table_ratio_mean.T.plot(kind='bar', ax=axes[2])
        axes[2].set_title('RT/ACC Ratio')
        axes[2].set_xlabel('tSOA Group')
        axes[2].set_ylabel('Ratio Value')
        axes[2].legend(title='Cue Type', loc='upper left', bbox_to_anchor=(1, 1))
        if len(subIDs) > 1:
            # 进行单侧配对t检验
            t_stat_rt, p_value_0 = ttest_rel(table_RT['Valid'], (table_RT['Short Invalid'] + table_RT['Long Invalid']) / 2, alternative='less')
            t_stat_acc, p_value_1 = ttest_rel(table_Acc['Valid'], (table_Acc['Short Invalid'] + table_Acc['Long Invalid']) / 2, alternative='greater')
            t_stat_ratio, p_value_2 = ttest_rel(table_ratio['Valid'], (table_ratio['Short Invalid'] + table_ratio['Long Invalid']) / 2, alternative='less')

            # 在相应子图标题下方显示p值
            axes[0].text(0.5, 0.95, f'p={p_value_0:.3f}', ha='center', va='bottom', color='red', fontsize=10, transform=axes[0].transAxes)
            axes[1].text(0.5, 0.95, f'p={p_value_1:.3f}', ha='center', va='bottom', color='red', fontsize=10, transform=axes[1].transAxes)
            axes[2].text(0.5, 0.95, f'p={p_value_2:.3f}', ha='center', va='bottom', color='red', fontsize=10, transform=axes[2].transAxes)

        # 设置整体的标题和布局
        plt.suptitle(f'Subjects: {subIDs}, Modality: {prefix_types}', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        # plt.tight_layout()
        plt.show()
            
# %% only correct RT
import pandas as pd
import glob
import matplotlib.pyplot as plt
import YLutilpy
import numpy as np

# 初始化设置
YLutilpy.default_img_set()
# for subIDs in [[2,3,5,7,8,9,10,11]]:  # 被试编号
for subIDs in [[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[2,3,4,5,6,7,8,9,10,11]]:  # 被试编号
# for subIDs in [[4]]:  # 被试编号
    for prefix_types in [['A'],['V']]:  # 组编号
        groupID = 1
        all_dfs = []

        # 加载数据
        # 初始化存储每个被试的平均值
        rt_means = {}
        acc_means = {}
        ratio_means = {}
        for prefix_type in prefix_types:  # 可扩展到其他模态如'A'
            for subID in subIDs:
                prefix = f"../Data/{prefix_type}_Result_G{groupID}_Sub{subID}_"
                files = glob.glob(prefix + "*.csv")
                if not files:
                    continue
                filename = files[0]
                df_tmp = pd.read_csv(filename)
                df_tmp = df_tmp[df_tmp['ID'] > 0]  # 保留有效ID
                df_tmp['subID'] = subID
                df_tmp['modality'] = prefix_type
                # 去除RT 大于或小于2.5个标准差的行
                df_tmp = df_tmp[(df_tmp['RT'] < df_tmp['RT'].mean() + 2.5 * df_tmp['RT'].std()) & 
                                (df_tmp['RT'] > df_tmp['RT'].mean() - 2.5 * df_tmp['RT'].std()) |
                                (df_tmp['RT'].isna())]  # 保留RT为NaN的行（错误反应）
                all_dfs.append(df_tmp)


        if not all_dfs:
            raise FileNotFoundError("No data files found.")

        df = pd.concat(all_dfs, ignore_index=True)

        df['tSOA_group'] = df.groupby('cueType')['tSOA'].rank(method='dense').astype(int)    
        # 将tSOA_group映射成字符标签
        tsoa_group_mapping = {1: 'Short Invalid', 2: 'Valid', 3: 'Long Invalid'}
        df['tSOA_group'] = df['tSOA_group'].map(tsoa_group_mapping)   

        # 将数据转换为表格形式
        table_RT = df[df['judge'] == 1].groupby(['subID', 'tSOA_group','cueType'])['RT'].mean().unstack(['tSOA_group','cueType'])
        table_Acc = df.groupby(['subID', 'tSOA_group','cueType'])['judge'].mean().unstack(['tSOA_group','cueType'])
        table_ratio = np.divide(table_RT,table_Acc)
        
        # 计算正确反应的RT（仅统计judge==1的RT）
        correct_df = df[df['judge'] == 1]  # 筛选正确反应
        rt_grouped = correct_df.groupby(['tSOA_group', 'cueType'])['RT'].mean().reset_index()

        # 计算正确率（所有试次）
        acc_grouped = df.groupby(['tSOA_group', 'cueType'])['judge'].mean().reset_index()

        # 创建透视表
        pivot_rt = rt_grouped.pivot(index='tSOA_group', columns='cueType', values='RT')
        pivot_acc = acc_grouped.pivot(index='tSOA_group', columns='cueType', values='judge')

        # 排序与颜色设置
        cue_order = ['AP2', 'PP', 'AP1', 'AU']
        color_dict = {
            'AP1': '#FFE699',
            'AP2': '#FFD966',
            'PP': '#42A62A',
            'AU': "#EC7025",
        }

        # 按顺序重排透视表
        pivot_rt = pivot_rt[cue_order]
        pivot_acc = pivot_acc[cue_order]
        pivot_ratio = pivot_rt / pivot_acc  # RT与正确率比值

        # 绘制图表
        fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=300)
        
        # 调整tSOA_group的顺序为 Short Invalid -> Valid -> Long Invalid
        tsoa_order = ['Short Invalid', 'Valid', 'Long Invalid']

        # 重新排序pivot表
        pivot_rt = pivot_rt.reindex(tsoa_order)
        pivot_acc = pivot_acc.reindex(tsoa_order)
        pivot_ratio = pivot_ratio.reindex(tsoa_order)
        # RT子图
        pivot_rt.plot(kind='bar', ax=axes[0], color=[color_dict[c] for c in cue_order])
        axes[0].set_title('Mean Reaction Time (Correct Trials Only)')
        axes[0].set_xlabel('tSOA Group')
        axes[0].set_ylabel('RT (s)')
        axes[0].get_legend().remove()

        # ACC子图
        pivot_acc.plot(kind='bar', ax=axes[1], color=[color_dict[c] for c in cue_order])
        axes[1].set_title('Accuracy (All Trials)')
        axes[1].set_xlabel('tSOA Group')
        axes[1].set_ylabel('Accuracy')
        axes[1].get_legend().remove()

        # 比值子图
        pivot_ratio.plot(kind='bar', ax=axes[2], color=[color_dict[c] for c in cue_order])
        axes[2].set_title('RT/ACC Ratio')
        axes[2].set_xlabel('tSOA Group')
        axes[2].set_ylabel('Ratio Value')
        axes[2].legend(title='Cue Type', loc='upper left', bbox_to_anchor=(1, 1))

        # 如果subIDs长度大于1，计算逐个被试的平均值并进行配对t检验
        if len(subIDs) > 1:
                # 进行配对t检验
            t_stat_rt, p_value_0 = ttest_rel(table_RT[('Valid', 'AU')], table_RT[('Valid', 'PP')])
            t_stat_acc, p_value_1 = ttest_rel(table_Acc[('Valid', 'AU')], table_Acc[('Valid', 'PP')])
            t_stat_ratio, p_value_2 = ttest_rel(table_ratio[('Valid', 'AU')], table_ratio[('Valid', 'PP')])    # 在柱状图上方标注P值
            
            for ax, p_value, label in zip(
                [axes[0], axes[1], axes[2]],
                [p_value_0, p_value_1, p_value_2],
                ['RT vs ACC', 'RT vs ACC', 'RT vs Ratio']
                ):
                ax.text(
                    0.5, 0.95, f'p={p_value:.3f}',
                    ha='center', va='bottom', color='red', fontsize=10,
                    transform=ax.transAxes
                )

        # 设置整体的标题和布局
        plt.suptitle(f'Subjects: {subIDs}, Modality: {prefix_types}', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        # plt.tight_layout()
        plt.show()

# %% 
import pandas as pd
import glob
import matplotlib.pyplot as plt
import YLutilpy
import numpy as np
from scipy.stats import ttest_ind
YLutilpy.default_img_set()
import matplotlib.pyplot as plt
YLutilpy.default_img_set()

subIDs = [1,2]  # 这里填你要分析的被试编号
groupID = 1

all_dfs = []
for subID in subIDs:
    for prefix_type in ['V','A']:
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

# %% 浓缩版分析
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

# %% threshold stage analysis
# 
import YLutilpy
YLutilpy.default_img_set()
import pandas as pd
import numpy as np
import glob 
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt
#取消网格线
plt.rcParams['axes.grid'] = False
subIDs = [1,2]
groupIDs = [1]
modalitys = ['A']  # 或 'V','A'，根据需要选择
# 构建文件名前缀
for subID in subIDs:
    for groupID in groupIDs:
        for modality in modalitys:
    
            prefix = f"../Data/{modality}_Result_G{groupID}_Sub{subID}_"
            files = glob.glob(prefix + "*.csv")
            if not files:
                continue
            filename = files[0]
            df_tmp = pd.read_csv(filename)
            aveACC = df_tmp[df_tmp['ID'] > 0]['judge'].mean()

            df_tmp = df_tmp[df_tmp['ID'] <= 0]  # 保留ID<=0且judge==1的行
            
            # 若存在tTilt列，删除tTilt等于0的行
            if 'tTilt' in df_tmp.columns:
                df_tmp = df_tmp.dropna(subset=['tTilt'])
                df_tmp = df_tmp[df_tmp['tTilt'] != 0]
            # 若存在tFreq列，删除tFreq等于0或者nan的行
            elif 'tFreq' in df_tmp.columns:
                df_tmp = df_tmp.dropna(subset=['tFreq'])
                df_tmp = df_tmp[df_tmp['tFreq'] != 0]    # 绘制tgAmp列的变化曲线
            # 删除tgAmp列为nan的行
            df_tmp = df_tmp.dropna(subset=['t0'])
                

            # 绘制tgAmp列的变化曲线
            plt.figure(figsize=(10, 5),dpi=300)
            plt.plot(df_tmp.index+1, df_tmp['tgAmp'], marker='o', linestyle='-', color='blue', label='tgAmp')
            
            # 标注转折点
            y = df_tmp['tgAmp'].values
            extremum_count = 0
            changeIdx = []
            for i in range(1, len(y) - 1):
                # 寻找i之前第一个不与y[i]相等的点
                prev_idx = i - 1
                while prev_idx >= 0 and y[prev_idx] == y[i]:
                    prev_idx -= 1
                
                # 如果所有前面的点都相等，则不是极值点
                if prev_idx < 0:
                    continue

                # 寻找i之后第一个不与y[i]相等的点
                next_idx = i + 1
                while next_idx < len(y) and y[next_idx] == y[i]:
                    next_idx += 1
                
                # 如果所有后面的点都相等，则不是极值点
                if next_idx >= len(y):
                    continue

                # 检查当前点i是否是平台的最后一个点
                is_last_of_plateau = (y[i] != y[i+1])

                # 判断是否为极值点
                is_extremum = (y[prev_idx] < y[i] and y[i] > y[next_idx]) or \
                              (y[prev_idx] > y[i] and y[i] < y[next_idx])

                if is_extremum and is_last_of_plateau:
                    extremum_count += 1
                    plt.scatter(df_tmp.index.values[i]+1, y[i], c='red', s=100, zorder=5)
                    plt.text(df_tmp.index.values[i]+1, y[i] + 0.01 * (max(y) - min(y)), str(extremum_count), 
                             ha='center', va='bottom', fontsize=12, color='black', weight='bold')
                    # 同时显示index值
                    plt.text(df_tmp.index.values[i]+1, y[i] - 0.02 * (max(y) - min(y)), str(df_tmp.index.values[i]+1), 
                             ha='center', va='top', fontsize=10, color='black')
                    # 记录该点index
                    changeIdx.append(df_tmp.index.values[i]+1)
            
            
            # 标记最后一个点的tgAmp
            plt.scatter(df_tmp.index.values[i+1]+1, df_tmp['tgAmp'].values[i+1], c='green', s=100, zorder=5)
            plt.text(df_tmp.index.values[i+1]+2, df_tmp['tgAmp'].values[i+1], f'Thr: {df_tmp["tgAmp"].values[i]:.5f}',
                     ha='left', va='center', fontsize=12, color='black', weight='bold')
            # 标记倒数第4个转折点至最末所有数据点的平均值
            
            last_four_avg = df_tmp['tgAmp'].values[df_tmp.index.values>=changeIdx[-4]].mean()
            
            
            
            # 绘制一条水平线横跨changeIdx[-4]到end
            start_x = changeIdx[-4]
            end_x = df_tmp.index.values[-1] + 1
            plt.hlines(y=last_four_avg, xmin=start_x, xmax=end_x, color="#DC6300", linestyle='--', label=f'Last 4 Avg: {last_four_avg:.5f}', linewidth=3.5)            
            
            plt.xlabel('Trial Index')
            plt.ylabel('tgAmp')
            plt.title(f'Threshold Stage - Sub{subID}, Modality: {modality}, Avg ACC: {aveACC:.3f}')            
            plt.legend()
            plt.show()
            # %% 
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

