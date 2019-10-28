from json import JSONDecodeError

import requests
from tqdm import tqdm
import pandas as pd
from lxml import html
import numpy as np
import os
from functools import lru_cache
import logging
logging.basicConfig(level=logging.INFO)
from file_cache.utils.util_log import timed, logger, timed_bolck
from file_cache.cache import file_cache
import time
from glob import glob

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


@lru_cache()
#@timed()
def get_mol_detail(smiles_id):
    session = get_session()
    link = f'https://ochem.eu/molecule/profile.do?depiction={smiles_id}&render-mode=popup'
    # print(link)
    try:
        response = session.get(link)
        if response.text.find('The molecule profile is unavailable') >=0 :
            logger.info(f'This is an empty molecule:{smiles_id}')
            return {'smiles_id':smiles_id, 'SMILES':'unavailable'}
    except Exception as e:
        logger.warning(e)
        return get_mol_detail(smiles_id)

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
        logger.info(f'Warning: can not get smiles from link: {link}')
        return {'smiles_id': smiles_id, 'SMILES': '0'}
    else:
        smiles_att['smiles_id'] = smiles_id
        return smiles_att


@timed()
def process_one_page(pagenum, property_id, pagesize):
    sleep_time = np.random.randint(1, 3)
    time.sleep(np.random.randint(sleep_time))

    args_local = locals()
    total_num, property_name = get_total_cnt(property_id)
    total_page = int(np.ceil(int(total_num) / int(pagesize)))


    fold = f'./output/ochem/{property_id:03}_{property_name}'
    fold = fold.replace(' ', '_')
    fold = fold.replace('(', '_')
    fold = fold.replace(')', '')
    os.makedirs(fold, exist_ok=True)
    df_file = f'{fold}/{property_id:03}_{pagesize}_{pagenum:04}.h5'

    if os.path.exists(df_file):
        exist_df = pd.read_hdf(df_file, 'mol')
        logger.info(f'Already have df:{exist_df.shape} in file:{df_file}')
        return exist_df
        # columns_cnt_base = 7
        # if len( exist_df.columns) > columns_cnt_base:
        #     logger.info(f'Already have df:{exist_df.shape} with col printableValue_ex  {property_name}/{property_id}')
        #     return exist_df
        # if len(exist_df.columns) == columns_cnt_base and not has_printvalue_ex(property_id):
        #     logger.info(f'Don not find printableValue_ex for {property_name}/{property_id} df:{exist_df.shape}')
        #     return exist_df
        # logger.warning(f'Rerun the batch process for {property_name}/{property_id} df:{exist_df.shape}, columns_cnt_base{columns_cnt_base}')

    res = get_request(property_id, pagenum, pagesize)

    if res.get('filters') is None or res.get('filters').get('filter') is None:
        if res is None:
            source = 1
        elif 'message' in res  \
                and  'message' in res.get('message') \
                and 'You are not logged in' in res.get('message').get('message'):
            raise Exception('Session is timeout')
        elif res.get('filters') is None:
            print(res)
            source = 2
        elif res.get('filters').get('filter') is None:
            source = 3
        logger.info(f'Error: Get none from process_one_page({source}):{args_local}')
        return process_one_page(pagenum, property_id, pagesize)


    #property_id = int(prop['value'])

    list_len = len(res.get('list').get('exp-property'))

    logger.info(f'There are {list_len}/{total_num} item in the res for page:{pagenum}/{total_page}, pagesize:{pagesize}')

    session = get_session()

    mol_list = []

    item_list = res.get('list').get('exp-property')

    for item in tqdm(item_list, f'{property_name}:{pagenum}/{total_page}'):
        try:
            step=1
            RecordID = item.get('id')
            step = 2
            printableValue = item.get('printableValue')
            step = 3
            doi = item.get('article').get('doi')
            step = 4
            #abb = item.get('article').get('journal').get('abbreviation')
            MoleculeID = item.get('molecule').get('mp2')
            step = 5
            # molWeight  =     item.get('molecule').get('molWeight')
            smiles_id = item.get('molecule').get('id')
            step = 6
            #smiles_att = get_mol_detail(smiles_id)
            step = 7

            right_att = {}
            if 'conditions' in item and 'property-value' in item.get('conditions'):
                prpperty_value = item.get('conditions').get('property-value')
                if isinstance(prpperty_value, list):
                    for one_property in prpperty_value :
                        key = one_property.get('property').get('name')
                        val = one_property.get('printableValue')
                        right_att[key] = val
                else:
                    key = prpperty_value.get('property').get('name')
                    val = prpperty_value.get('printableValue')
                    right_att[key] = val


        except Exception as e:
            logger.exception(e)
            logger.error(f'Parse error for :{RecordID}, step:{step}, {(pagenum, property_id, pagesize)}')
            logger.error(item)
            return e

        mol = {
            'RecordID': RecordID,
            'printableValue': printableValue,
            'doi': doi,
            #'abb': abb,
            'MoleculeID' : MoleculeID,
            # 'molWeight':molWeight,
            'smiles_id' : smiles_id,
            'property_name' : property_name,
            'property_id' : property_id,
            #'printableValue_ex' : printableValue_ex,

        }
        mol.update(right_att)

        #mol.update(smiles_att)
        mol_list.append(mol)

    df = pd.DataFrame(mol_list)

    columns = ['RecordID', 'MoleculeID', 'SMILES', 'doi',
               #'abb',
               'printableValue', 'property_name', 'property_id',
               'Name', 'Formula', 'InChI Key', 'Molecular Weight', 'smiles_id']

    for key in mol.keys():
        if key not in columns:
            columns.append(key)



    columns = [col for col in columns if col in df]
    df[columns].to_hdf( df_file, 'mol')
    logger.info(f'file save to {df_file}')

    return df


@timed()
@lru_cache()
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
    r = session.post(url=url, data=data, headers=headers1, timeout=200)

    # get the request info:
    # List size, and list property_name, property_id
    import json
    try:
        res = json.loads(r.text.encode('utf-8'))
    except JSONDecodeError as e:
        logger.exception(e)
        logger.info(f'Error for paras:{property_id, pagenum, pagesize}')
        logger.info(res[:100])
        raise e

    return res

@lru_cache()
def get_total_cnt(property_id):
    res = get_request(property_id, 1, pagesize=5)
    total_num = res.get('list').get('size')

    prop = [item for item in res.get('filters').get('filter') if item.get('name') == 'property']

    prop = dict(prop[0])

    property_name = prop['title']
    return total_num, property_name

@lru_cache()
def has_printvalue_ex(property_id):
    res = get_request(property_id, 1, pagesize=5)
    item_list = res.get('list').get('exp-property')
    return 'property-value' in item_list[0].get('conditions')

@timed()
def process_one_item(property_id, thread_num=3):
        pagesize = 50
        total_num, property_name = get_total_cnt(property_id)
        total_page = int(np.ceil(int(total_num) / int(pagesize)))

        # for pagenum in range(1, total_page+1):
        #     process_one_page(property_id, pagenum=pagenum, pagesize=pagesize)

        #from multiprocessing import Pool as ThreadPool  # 进程
        from multiprocessing.dummy import Pool as ThreadPool  # 线程
        pool = ThreadPool(thread_num)

        from functools import partial
        process_one_page_ex = partial(process_one_page,  property_id=property_id, pagesize=pagesize)

        df_list = pool.map(process_one_page_ex, range(1, total_page+1), chunksize=1)

        return pd.concat(df_list)




#@timed()
def fill_smiles(property_id, thread_num=3):
    res = get_request(property_id, 1, 5)

    prop = [item for item in res.get('filters').get('filter') if item.get('name') == 'property']

    prop = dict(prop[0])

    property_name = prop['title']
    property_name = property_name.replace(' ', '_')
    property_name = property_name.replace('(', '_')
    property_name = property_name.replace(')', '')

    property_id_query = int(prop['value'])

    if property_id != property_id_query :
        return f'{property_id} vs {property_id_query}'


    reg = f'./output/ochem/{property_id:03}_{property_name}/*.h5'

    file_list = glob(reg)
    file_list = sorted(set(file_list))

    logger.info(f'find {len(file_list)} file with :{reg}')
    from multiprocessing.dummy import Pool as ThreadPool  # 线程
    pool = ThreadPool(thread_num)

    df_list = pool.map(fill_smiles_file, file_list, chunksize=1)
    return pd.concat(df_list)


@timed()
def fill_smiles_file(file):
    file_name = os.path.basename(file)
    try:
        try:
            print(1)
            df = pd.read_hdf(file, 'mol').drop_duplicates('smiles_id')
            logger.warning(f'No smiles in file:{file_name}')
            print(2)
            with pd.HDFStore(file, mode='r') as store:
                keys = store.keys()
                if '/smiles' in keys:
                    exist_smiles = pd.read_hdf(file, 'smiles')
                else:
                    exist_smiles = pd.DataFrame(columns=['smiles_id'])
        except Exception as e :
            logger.exception(e)
        all_smiles = get_smiles_all()
        exist_smiles = pd.concat([exist_smiles, all_smiles])

        exist_smiles = exist_smiles.loc[exist_smiles.SMILES.fillna('').str.len() > 1]

        todo = df.loc[~df.smiles_id.isin(exist_smiles.smiles_id)]

        if len(todo) > 0:
            tqdm.pandas(desc=file_name)
            smiles_list = todo.smiles_id.progress_apply(get_mol_detail).to_list()
            smiles_df_new = pd.DataFrame(smiles_list)
            smiles_df_new = smiles_df_new.dropna(how='all')

            smiles_df = pd.concat([exist_smiles, smiles_df_new])

            #Pick up simles for current DF
            smiles_df = smiles_df.loc[smiles_df.smiles_id.astype(str).isin(df.smiles_id.astype(str))]

            smiles_df.to_hdf(file, 'smiles')

            logger.info(f'Find {len(todo)} smiles from {len(smiles_df_new)} todo, exist smiles:{len(exist_smiles)}, total mol:{len(df)}')
            logger.info(f'Current status:{len(smiles_df)}/{len(df)}, {file}')


            if len(smiles_df) == 0:
                logger.warning(f'Can not find similes for file:{file}')
        else:
            print('='*100)
            logger.info((type(df), type(exist_smiles)))
            logger.info((df.dtypes, exist_smiles.dtypes))
            print(type(exist_smiles.smiles_id), type(exist_smiles['smiles_id']), exist_smiles.head().smiles_id, df.smiles_id.to_list())
            smiles_df = exist_smiles.loc[exist_smiles['smiles_id'].astype(str).isin(df['smiles_id'].astype(str))]
            logger.warning(f'No similes need to pull for file:{file}, todo:{len(todo)}, cur_similes:{len(smiles_df)}, all_smiles:{len(all_smiles)},exist_smiles:{len(exist_smiles)}')
            #Pick up simles for current DF
        return smiles_df
    except Exception as e:
        logger.exception(e)
        logger.error(f'There is error when process:{file}')
        return pd.DataFrame()

@timed()
@file_cache()
def get_smiles_all():
    df_list = []
    df_list.append(pd.DataFrame())

    file_list = list(glob('./output/ochem/*/*.h5'))
    for file in file_list:
        try:
            smiles = pd.read_hdf(file, 'smiles')
            df_list.append(smiles)
            #logger.info(f'Find {len(smiles)} smiles in file:{file}')
        except Exception as e :
            #logger.exception(e)
            #logger.warning(f'No smiles df in file:{file}')
            pass
    df = pd.concat(df_list)

    df = df.loc[df.SMILES.fillna('').str.len()>1]
    print(df.shape, df.dtypes, type(df.smiles_id), df.head().smiles_id)
    return df.drop_duplicates(['smiles_id'])

@timed()
@file_cache(type='h5')
def get_feature(fold='./output/ochem/156_Plasma_protein_binding'):
    file_list = sorted(list(glob(f'{fold}/*.h5')))

    mol_list = []
    for file in tqdm(file_list, fold):
        mol = pd.read_hdf(file, 'mol')
        smiles = get_smiles_all()

        mol = mol.merge(smiles, how='left', on='smiles_id')
        mol_list.append(mol)
    return pd.concat(mol_list).df.drop_duplicates()

def get_feature_final():
    from glob import glob
    os.makedirs('./final/', exist_ok=True)
    for fold in list(glob('./output/*/*')):
        try:
            df = get_feature(fold)
            file = f'final/{os.path.basename(fold)}_{df.shape[0]}_{df.shape[1]}.h5'
            df.to_hdf(file, 'key')
            print(file, df.shape)
        except Exception as e:
            print(f'Error:{fold} ')


if __name__ == '__main__':
    #218 CYP450_modulation
    import fire
    fire.Fire()
    #fill_smiles_file('./output/ochem/004_LogL_water/004_50_0001.h5')
    #process_one_item(218, thread_num=3)

    #get_smiles_all()