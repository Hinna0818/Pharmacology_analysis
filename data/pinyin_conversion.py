import pyreadr
from pypinyin import lazy_pinyin
#print(''.join(lazy_pinyin('地龙')))

## 读取数据
data = pyreadr.read_r('./input your herb data.rda')
result = data['your herb data name']

# 将非字符串类型的数据转换为字符串
result['Herb_cn_name'] = result['Herb_cn_name'].astype(str)

# 将中文名转换为拼音名
result['Herb_pinyin_name'] = result['Herb_cn_name'].apply(lambda x: ''.join(lazy_pinyin(x)))

pyreadr.write_rdata('./ouput_data.rda', result)


