# dig

## Structure


- *dig*, python package folder
- *notebooks*, python notebook folder to demonstrate how can we use dig

## Install

Prequisite: [MDSplus](https://mdsplus.org/index.php?title=Downloads&open=1769734509572583850001&page=Software%2FDownloads)

```bash
# for clients
pip install dig-mhd 

# for developers
cd /to-somewhere-you-want-to-clone-the-repo
git clone https://github.com/WenyinWei/dig.git
cd dig # (not dig/dig)
pip install -e . # "-e" means editable
```

## Use Python in Matlab

### Common Data Format

_*.npy_, _*.npz_ file formats are common in python/numpy community; while _*.mat_ file format is common in matlab community. Matlab has professional [documents](https://ww2.mathworks.cn/products/matlab/matlab-and-python.html) to demonstrate how to use Python from Matlab and vice.

```python
from scipy.io import loadmat
annots = loadmat('*.mat')
```

```matlab
% In Matlab, npy-matlab module could help you read and write npy format files, https://github.com/kwikteam/npy-matlab
a = rand(5,4,3);
writeNPY(a, 'a.npy');
b = readNPY('a.npy');
```

```matlab
% In Matlab, configure Python environment https://ww2.mathworks.cn/help/matlab/ref/pyenv.html
pyenv % 返回当前默认的 Python 环境
pyenv('Version', 3.8) % 寻找当前电脑上的对应版本的 Python 环境
pyenv('Version', path-to-your-python-executable) % 寻找当前电脑上的对应路径的 Python 环境

if count(py.sys.path,'path-to-your-package') == 0
    insert(py.sys.path,int32(0),'path-to-your-package');
end

py.importlib.import_module('dig')
py.importlib.reload('dig')
returned_tuple = py.dig.utilvar.get_EAST_EFIT_BR_BZ_Bt(shotnum, py.list({tpoint1, tpoint2, ...}));
R = returned_tuple{1};
Z = returned_tuple{2};
tpoints = returned_tuple{3};
BRs = returned_tuple{4};
BZs = returned_tuple{5};
Bts = returned_tuple{6};
```

```python
# Rewrite all npy files in mat format, https://ww2.mathworks.cn/matlabcentral/answers/444998-how-do-read-npy-files-in-matlab#answer_624742
from scipy.io import savemat
import numpy as np
import glob
import os

npzFiles = glob.glob("*.npz")
for f in npzFiles:
    fm = os.path.splitext(f)[0]+'.mat'
    d = np.load(f)
    savemat(fm, d)
    print('generated ', fm, 'from', f)
```