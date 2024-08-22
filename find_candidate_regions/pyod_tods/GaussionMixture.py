import numpy as np
from sklearn.mixture import GaussianMixture

##
# 加载数据
file_path = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/tods/test_new.bed"
raw_data = np.loadtxt(file_path)
data = np.expand_dims(raw_data[:586177], axis=1)
# 选择高斯分布的个数,这里取3个
n_components = 3

# 构建GMM模型,covariance_type选择full更加灵活
gmm = GaussianMixture(n_components=n_components, covariance_type='full', means_init=[[3], [30], [70]], max_iter=200)

# 训练模型
gmm.fit(data)

# 预测类别标签
labels = gmm.predict(data)

# 获取高斯分布的参数:均值、方差、权重 
means = gmm.means_ 
covs = gmm.covariances_
weights = gmm.weights_  

ls = []
# 判断区域特征
for i in range(n_components): 
    print("component {}".format(i))    
    # 遍历每个高斯分布
    idx = np.where(labels==i)[0]   # 当前高斯对应的所有数据点下标
    
    region = []  
    start = -1
    end = -1
    for j in idx:
        if start > 0:
            if j > end + 1:
                region.append([start, end])
                start = j
                end = j
            else:
                end = j
        else:
            start = j
            end = j
            region = []
    for reg in region:
        print(reg, end=",")

    # 这些数据点构成一个区域  
    # 通过高斯的参数判断区域特征    
    #   if means[i] > 3 and covs[i] < 1:  
    if weights[i] > 0.8:    # 权重最大的部分
        print('特征区域:',idx)  # 均值较高,方差较小的区域
        print(weights[i], means[i], covs[i])
    else: 
        print('噪声区域:',idx) # 权重较小的噪声区域
        print(weights[i], means[i], covs[i])
        ls.extend(list(idx))

# 我们还可以通过后验概率进行软聚类并判断区域特征  
probs = gmm.predict_proba(data) 
for i in range(n_components):
    idx = np.where(probs[:,i] > 0.5)[0] # 后验概率大于0.5的点
    # 通过probs[:,i]的均值判断该区域的特征 

region = []
start = -1
end = -1
for j in ls:
    if start > 0:
        if j - end > 20 :
            region.append([start, end])
            start = j
            end = j
        else:
            end = j
    else:
        start = j
        end = j
        region = []
f = open("/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/tods/out.bed", "w")
for reg in region:
    f.write("{}\t{}\t{}\n".format("NC_000019.10", reg[0] * 100, reg[1] * 100))
    # print(reg, end=",")
f.close()