import pandas as pd
import numpy as np
import pubchempy as pbp
import re
import tqdm
import threading
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import BoundaryNorm


def get_compounds(indices, result='cid', identifier='name') -> list: 
    data = {i:j for i,j in zip(indices, range(len(indices)))}
    threads = []

    def get_compound(index):
        while True:
            if identifier == 'name':
                name = index
                try:
                    compound = pbp.get_compounds(name, 'name')[0]
                    break
                except IndexError:
                    try:
                        compound = pbp.get_compounds(name.split(' ')[0], 'name')[0]
                        break
                    except IndexError:
                        try:
                            compound = pbp.get_compounds(name.split('（')[0], 'name')[0]
                            break
                        except IndexError:
                            print('{} is not found'.format(name))
                            return None
                        except pbp.PubChemHTTPError:
                            time.sleep(0.1)
                    except pbp.PubChemHTTPError:
                        time.sleep(0.1) 
                except pbp.PubChemHTTPError:
                    time.sleep(0.1) 
            else:
                try:
                    compound = pbp.get_compounds(index, identifier)[0]
                    break
                except IndexError:
                    print('this compound is not find, please check out your identifier.')
                except pbp.PubChemHTTPError:
                    time.sleep(0.1)
        if result == 'cid':
            data[index] = compound.cid
        elif result == 'fingerprint':
            data[index] = compound.fingerprint

    for i,s in enumerate(indices):
        threads.append(threading.Thread(target=get_compound, args=(s,)))
        threads[i].start()
    for t in threads:
        t.join()
    return list(data.values())

def notxa0(names):
    res = []
    for s in names:
        res.append(s.replace("\xa0", " "))
    return res

def add_data(df, datadict):
    cids = datadict['cid']
    for i in range(len(cids)):
        if cids[i] in df['cid']:
            j = df.loc[cids[i], 'cid']
            for s in datadict.keys():
                df[s][j] = datadict[s][i]

def strtonp(s):
    numbers = re.findall('[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', s)
    numbers = np.asarray([float(i) for i in numbers])
    return numbers

def get_CI(df, cell_lines=None):
    CI = []
    titles= df.columns.to_list()[4:-1]
    if cell_lines == None:
        for i in titles:
            for j in df[i]:
                if not pd.isna(j):
                    CI.append(strtonp(j))
    elif cell_lines in titles:
        for i in df[cell_lines]:
            if not pd.isna(i):
                    CI.append(strtonp(i))
    else:
        print('{} is wrong'.format(cell_lines))
    return CI

def divide(x,y):
    data = np.array(list(zip(x,y)))
    p1, p2, p3, p4 = [], [], [], []
    for i in data:
        if i[0] <= 1 and i[1] <= 1:
            p1.append(i)
        elif i[0] >1 and i[1] > 1:
            p3.append(i)
        elif i[0] > 1:
            p2.append(i)
        else:
            p4.append(i)
    return p1, p2, p3, p4

def trans_circle(array, center=(0, 0)):
    array = np.asarray(array)
    l = np.sqrt((array[:, 0]-center[0])**2+(array[:,1]-center[1])**2)
    j = np.arctan2(array[:, 1]-center[0], array[:, 0]-center[0])
    return j, l

def drawing_3points(CI, mode=0, name=None, theindex=None):
    x, y, z = zip(*CI)
    line = np.log(np.linspace(min(min(y), min(z)), max(max(y), max(z)), 5))
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.dpi = 1024
    if mode == 0:
        formatter = FuncFormatter(lambda x, pos: '{:.2f}'.format(np.exp(np.exp(x)-1)-1),)
        ax.yaxis.set_major_formatter(formatter)
        p1, p2, p3, p4 = divide(y, z)
        x_, y_ = trans_circle(list(zip(y, z)), (1, 1))
        boundaries = list(np.linspace(min(x), 0.99, 99)) + list(np.linspace(0.99, 1.01, 101)) + list(np.linspace(1.01, max(x), 100))
        scatter = ax.scatter(x_, np.log(np.log(y_+1)+1), s=2, c = x, cmap='viridis', norm=BoundaryNorm(boundaries, ncolors=300, clip=True))
        cbar = fig.colorbar(scatter, ax=ax, ticks=[])
        cbar.set_ticks([min(x), 1, max(x)])
        # cbar.set_ticklabels(['synergy', 'additivity', 'antagonism'])
        y_ = np.log(y_+1)
        cbar.set_label('0.1uM NCK_CI')  # 设置颜色条的标签
        ax.plot([0]*5, np.linspace(0, max(np.log(y_+1)), 5), '--', lw=0.5, color='#1f77b4')
        ax.plot([np.pi/2]*5, np.linspace(0, max(np.log(y_+1)), 5), '--', lw=0.5, color='#1f77b4')
        ax.plot([-np.pi/2]*5, np.linspace(0, max(np.log(y_+1)), 5), '--', lw=0.5, color='#1f77b4')
        ax.text(-2.95, max(np.log(y_+1))-0.1, '$CI_1<1,CI_2<1$\nnumbers:{}'.format(len(p1)),size=8)
        ax.text(-2.95+np.pi/2, max(np.log(y_+1))-0.1, '$CI_1>1,CI_2<1$\nnumbers:{}'.format(len(p2)),size=8)
        ax.text(-2.95+np.pi, max(np.log(y_+1))-0.1, '$CI_1>1,CI_2>1$\nnumbers:{}'.format(len(p3)),size=8)
        ax.text(-2.95+3*np.pi/2, max(np.log(y_+1))-0.1, '$CI_1<1,CI_2>1$\nnumbers:{}'.format(len(p4)),size=8)
        ax.set_xlim(-3.15, 3.15)
        ax.set_title(''.format(name))
        ax.set_xlabel('angle')
        ax.set_ylabel('distance')
    elif mode==1:
        formatter = FuncFormatter(lambda x, pos: '{:.2f}'.format(np.exp(x)),)
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_formatter(formatter)
        p1, p2, p3, p4 = divide(y, z)
        x_, y_ = np.log(y), np.log(z)
        boundaries = list(np.linspace(min(x), 0.99, 99)) + list(np.linspace(0.99, 1.01, 101)) + list(np.linspace(1.01, max(x), 100))
        scatter = ax.scatter(x_, y_, s=2, c = x, cmap='viridis', norm=BoundaryNorm(boundaries, ncolors=300, clip=True))
        cbar = fig.colorbar(scatter, ax=ax, ticks=[])
        cbar.set_ticks([min(x), 1, max(x)])
        # cbar.set_ticklabels(['synergy', 'additivity', 'antagonism'])
        cbar.set_label('0.1uM NCK_CI')  # 设置颜色条的标签
        ax.plot([0]*5, line, '--', lw=0.5, color='#1f77b4')
        ax.plot(line, [0]*5, '--', lw=0.5, color='#1f77b4')
        ax.plot([min(line), max(line)], [min(line), max(line)], '--', lw=0.5, color='red')
        ax.set_title(''.format(name))
        ax.set_xlabel('1uM NCK')
        ax.set_ylabel('10uM NCK')
    if theindex != None:
        plt.savefig('../组会/{}.png'.format(theindex), dpi=1024)
        plt.close()



if __name__ == "__main__":
    import re

    # 定义CAS号的正则表达式
    cas_regex = r'\b[0-9]-[0-9A-Z]{2,}\b'

    # 测试字符串
    test_string = "这个化合物的CAS号是123-45-6，另一个化合物的CAS号是12345-67-8，还有一个是123456-78-9。"

    # 使用正则表达式查找所有匹配的CAS号
    matches = re.findall(cas_regex, test_string)

    # 打印匹配结果
    print(matches)