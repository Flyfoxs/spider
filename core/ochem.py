import requests
from tqdm import tqdm
import pandas as pd
from lxml import html
import numpy as np
import os
from functools import lru_cache
import logging
logging.basicConfig(level=logging.INFO)
from file_cache.utils.util_log import timed
from file_cache.cache import file_cache


@lru_cache()
def get_session():
    username='mindrank'
    password='mindrank1239'

    # 通过Session类新建一个会话
    session = requests.Session()
    post_url = 'https://ochem.eu/login/login.do'
    # 往下使用requests的地方，直接使用session即可，session就会保存服务器发送过来的cookie信息

    data = {
        'login': username, # 账号
        #'render-mode':'redirect',
        'pwd': password, # 密码

    }

    r = session.post(url=post_url, data=data)

    if r.text.find(username):
        return session
    else:
        raise Exception('login failed')

# function to get similes and other properties



def get_mol_detail(smiles_id):
    session = get_session()
    link = f'https://ochem.eu/molecule/profile.do?depiction={smiles_id}&render-mode=popup'
    # print(link)
    response = session.get(link)

    # print(response.text)
    tree = html.fromstring(response.text)
    # print(tree.get)
    # print(tree.text)
    tr_list = tree.xpath('//td/table[@class="properties"]/tr')
    # tr_list = tree.xpath('//td')

    # print(len(tr_list), tr_list[0].tag,  tr_list[0].keys(), dir(tr_list[0]) )
    # print(tr_list[2].text_content())
    smiles_att = {}
    for item in tr_list:
        # print('='*100)
        # print(item.text_content())
        td_list = item.xpath('td')

        key = td_list[0].text_content().strip()
        value = td_list[-1].text_content().strip()

        smiles_att[key] = value

    if len(smiles_att) == 0:
        print(f'Warning: can not get smiles from link: {link}')
    return smiles_att


@timed()
@file_cache()
def process_one_page(pagenum, property_id, pagesize):
    total_num = get_total_cnt(property_id)
    total_page = int(np.ceil(int(total_num) / int(pagesize)))

    res = get_request(property_id, pagenum, pagesize)
    prop = [item for item in res.get('filters').get('filter') if item.get('name') == 'property']
    prop = dict(prop[0])

    property_name = prop['title']
    property_id = int(prop['value'])

    list_len = len(res.get('list').get('exp-property'))

    print(f'There are {list_len}/{total_num} item in the res for page:{pagenum}/{total_page}, pagesize:{pagesize}')

    session = get_session()

    mol_list = []

    item_list = res.get('list').get('exp-property')

    for item in tqdm(item_list, f'{property_name}:{pagenum}/{total_page}'):
        RecordID = item.get('id')
        printableValue = item.get('printableValue')

        doi = item.get('article').get('doi')
        abb = item.get('article').get('journal').get('abbreviation')

        MoleculeID = item.get('molecule').get('mp2')
        # molWeight  =     item.get('molecule').get('molWeight')
        smiles_id = item.get('molecule').get('id')

        smiles_att = get_mol_detail(smiles_id)

        mol = {
            'RecordID': RecordID,
            'printableValue': printableValue,
            'doi': doi,
            'abb': abb,
            'MoleculeID': MoleculeID,
            # 'molWeight':molWeight,
            'smiles_id': smiles_id,
            'property_name': property_name,
            'property_id':property_name,

        }

        mol.update(smiles_att)
        mol_list.append(mol)

    df = pd.DataFrame(mol_list)

    columns = ['RecordID', 'MoleculeID', 'SMILES', 'doi',
               'abb', 'printableValue', 'property_name', 'property_id',
               'Name', 'Formula', 'InChI Key', 'Molecular Weight', 'smiles_id']

    fold = f'./output/ochem/{property_name}'
    fold = fold.replace(' ', '_')
    df_file = f'{fold}/{property_id:03}_{pagesize}_{pagenum:04}.h5'

    os.makedirs(fold, exist_ok=True)
    df[columns].to_hdf( df_file, 'mol')
    print(f'file save to {df_file}')

    return df


@timed()
def get_request(property_id, pagenum, pagesize = 500):
    url = 'https://ochem.eu/epbrowser/list.do'
    headers1 = {
        'Host': 'ochem.eu',
        'Connection': 'keep-alive',
        'Accept': 'application/json, text/javascript, */*; q=0.01',
        'Origin': 'https://ochem.eu',
        'X-Requested-With': 'XMLHttpRequest',
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/77.0.3865.90 Safari/537.36',
        'Sec-Fetch-Mode': 'cors',
        'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8',
        'Sec-Fetch-Site': 'same-origin',
        'Referer': 'https://ochem.eu/epbrowser/show.do?property=1&render-mode=popup',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7',
        # 'Cookie': 'fontsize=100; JSESSIONID=5A0F8AE2892C665E72EA2401F74E0093',
    }

    session = get_session()
    #for pagenum in range(1):
    data = {
        'property': property_id,
        'approval-status': 'all',
        'basket-select': '-1',
        'operation-id': '78438932',
        'out': 'json',
        'pagenum': pagenum,
        'pagesize': pagesize,
        'start-operation': '1',
        'structure-basket': '-1',
        'structure-basket-dofilter': '0',
        'visibility': 'all',
        'xemistry-similarity-cutoff': '85',
        'xemistry-sstype': 'substructure',
    }
    r = session.post(url=url, data=data, headers=headers1)

    # get the request info:
    # List size, and list property_name, property_id
    import json
    res = json.loads(r.text)

    return res

@lru_cache()
def get_total_cnt(property_id):
    res = get_request(property_id, 1, pagesize=5)
    total_num = res.get('list').get('size')
    return total_num


@timed()
def process_one_item(property_id):
        pagesize = 500
        total_num = get_total_cnt(property_id)
        total_page = int(np.ceil(int(total_num) / int(pagesize)))

        # for pagenum in range(1, total_page+1):
        #     process_one_page(property_id, pagenum=pagenum, pagesize=pagesize)

        from multiprocessing import Pool as ThreadPool  # 进程
        from multiprocessing.dummy import Pool as ThreadPool  # 线程
        pool = ThreadPool(3)

        from functools import partial
        process_one_page_ex = partial(process_one_page,  property_id=property_id, pagesize=pagesize)

        df_list = pool.map(process_one_page_ex, range(1, total_page+1), chunksize=1)





if __name__ == '__main__':
    process_one_item(1)