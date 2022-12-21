import sys
sys.path.insert(1, '')
from boimmg.boimmgpy.read_model.case_study import read_treat_model
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




model=read_treat_model()



results=model[0]
class_count={'Phosphatidylcholine': 113, 'CDP-1,2-diacyl-sn-glycerol': 20, 'Phosphatidylethanolamine': 72, 'Phosphatidylinositol': 12, 'Phosphatidylglycerophosphate': 2, 'Phosphatidylglycerol': 5}
check_annotation=model[2]
annotations_before_implementation=list(check_annotation.values())
print(len(annotations_before_implementation))
print(annotations_before_implementation.count(True))
print((annotations_before_implementation.count(True)/len(annotations_before_implementation))*100)

for l in check_annotation.keys():
    if l in results.keys():
        check_annotation[l]=True

annotations_after_implementation=list(check_annotation.values())
print(annotations_after_implementation.count(True))
print((annotations_after_implementation.count(True)/len(annotations_after_implementation))*100)

# desvio padrao e media


listr = []

     
for value in results.values():
    listr.append(value)

for count, value in enumerate(listr):
        listr[count]=len(value)

one_hit=0
sum_hits=sum(listr)
for value in listr:
    if value==1:
        count+=1
   
data=pd.Series(listr)
print(data.describe())

# grafico
names = ['Annoted before implementation', 'Annoted after implementation']
values = [(annotations_before_implementation.count(True)/len(annotations_before_implementation))*100 , (annotations_after_implementation.count(True)/len(annotations_after_implementation))*100]

y_pos = np.arange(len(names))
plt.bar(y_pos, values)

# Create names on the x-axis
plt.xticks(y_pos, names)
plt.yticks(np.arange(0, 101, 10))
# Show graphic
plt.show()




listr=[]

num_hits=0
total_one_it=0
for value in results.values():
    listr.append(value)

for count,value in enumerate(listr):
    num_hits+=len(value)
    if len(value)==1:
        total_one_it+=1

print(num_hits,total_one_it)


print(class_count)
