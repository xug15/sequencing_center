# 童明汉 Results:

> project_190620_TMH-XNN #(第三批)

## 实验设计

|实验组|对照组|
|:-|:-|
|ck0|WT|
|3号 ko 小鼠 |3号 WT 小鼠 |
|15号 ko 小鼠 |15号 WT 小鼠 |

每个都 total sequencing, Ribo-seq 数据。

## Raw data process pipeline (Ding Jiyu Finished)
**Qquality**
![FC](../image/190620-XNN/afterqality.png)  

## rRNA 
* rRNA 含量都在30-45%之间，此项结果合格。

| Iterm       | RC15              | RC3              | RW15              | RW3               | TC15              | TC3               | TW15              | TW3               | 
|-------------|-------------------|------------------|-------------------|-------------------|-------------------|-------------------|-------------------|-------------------| 
| Total reads | 29117568          | 15945845         | 29117568          | 26538239          | 21647172          | 21085495          | 20839854          | 21286104          | 
| rRNA        | 12868091 (44.19%) | 6231774 (39.08%) | 12868091 (44.19%) | 11132644 (41.95%) | 5195548 (24.00%)  | 7424015 (35.21%)  | 6839769 (32.82%)  | 6633799 (31.16%)  | 
| non rRNA    | 16249477 (55.81%) | 9714071 (60.92%) | 16249477 (55.81%) | 15405595 (58.05%) | 16451624 (76.00%) | 13661480 (64.79%) | 14000085 (67.18%) | 14652305 (68.84%) | 

## Read length distribution

* Ribo-seq library length 分布以30nt为均值的正态分布，结果合格。对照的转录组建库长度分布主要分布在30nt,结果合格。

RC15 
![FC](./result_190620_XNN/RC15_length.png)  
RC3
![FC](./result_190620_XNN/RC3_length.png)  
RW15 
![FC](./result_190620_XNN/RW15_length.png) 
RW3 
![FC](./result_190620_XNN/RW3_length.png)  
TC15
![FC](./result_190620_XNN/TC15_length.png)  
TC3
![FC](./result_190620_XNN/TC3_length.png)  
TW15
![FC](./result_190620_XNN/TW15_length.png)  
TW3
![FC](./result_190620_XNN/TW3_length.png)  

## read to RNA DNA and Intron (Using Tophat and readsNumCal_intron_v3)

* read 大部分都比对到RNA中，结果合理。

| Iterm                                | RC15     | RC3      | RW15     | RW3      | TC15     | TC3      | TW15     | TW3      | 
|--------------------------------------|----------|----------|----------|----------|----------|----------|----------|----------| 
| unique mapped reads of RNA           |  2781456 |  2052673 |  3081580 |  3376747 |  3404711 |  3176944 |  3445547 |  3377738 | 
| unique mapped reads of Intron        |  1459873 |  1136971 |  1735785 |  1886293 |  477659  |  569024  |  359882  |  522544  | 
| unique mapped ambiguous reads of RNA |  327081  |  196881  |  492936  |  564206  |  299193  |  183668  |  308807  |  307879  | 
| unique mapped reads of DNA           |  612731  |  381812  |  679016  |  536901  |  3750711 |  2842941 |  3290770 |  3258154 | 


## Differentail translatons (Xtail) Results: (Xu Gang Finished)

### File list

* [xnnFC.pdf](./result_190620_XNN/xnnFC.pdf)  
* [xnnRs.pdf](./result_190620_XNN/xnnRs.pdf)  
* [xnn_results.txt](./result_190620_XNN/xnn_results.txt)  
* [xnnfc_results.txt](./result_190620_XNN/xnnfc_results.txt)  
* [xnnmerge.counter](./result_190620_XNN/xnnmerge.counter)  
* [xnnrs_results.txt](./result_190620_XNN/xnnrs_results.txt)  
* [xnnvolcano.pdf](./result_190620_XNN/xnnvolcano.pdf)  


[xnn_results.txt](../image/190620-XNN/xnn_results.txt)  

![FC](../image/190620-XNN/FC.png)  
[FC](../image/190620-XNN/xnnFC.pdf)  

![Rs](../image/190620-XNN/Rs.png)  
[Rs](../image/190620-XNN/xnnRs.pdf)  


![Volcano](../image/190620-XNN/volcano.png)  
[Volcano](../image/190620-XNN/xnnvolcano.pdf)  








