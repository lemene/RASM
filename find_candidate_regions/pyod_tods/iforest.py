
####
import pandas as pd

from tods import schemas as schemas_utils
from tods import generate_dataset, evaluate_pipeline

import numpy as np 
from tods.detection_algorithm import PyodIsolationForest

table_path = '/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/depths/test/test.bed'
# target_index = 6 # what column is the target
# metric = 'F1_MACRO' # F1 on both label 0 and 1

# Read data and generate dataset
df = pd.read_csv(table_path, sep="\t")
dataset = generate_dataset(df)

# Load the default pipeline
method = PyodIsolationForest(n_estimators=100, max_features=1.0, contamination=0.1, verbose=1)
# 训练模型
method.fit(dataset)

# 获得异常分数
scores = method.decision_function(dataset)

# 设置阈值并获得预测标签
labels = scores < 0   
# labels中1表示正常,-1表示异常
# 评估效果
# n_errors = (labels != true_labels).sum()
# accuracy = 100*(1 - n_errors/len(labels))  

# 打印检测结果
for i in range(len(scores)):
    if labels[i] == 1:
        print(f"正常点,坐标:{dataset[i,0]},得分{scores[i]}")
    else:
        print(f"异常点,坐标:{dataset[i,0]},得分{scores[i]}") 



# pipeline = 

# # Run the pipeline
# pipeline_result = evaluate_pipeline(dataset, pipeline)
# print(pipeline_result)