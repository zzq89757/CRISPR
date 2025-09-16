import mysql.connector
import csv
import re
import zipfile
import os

def clean_column_name(name: str) -> str:
    """ 清洗列名，替换特殊字符为下划线 """
    return re.sub(r'\W|^(?=\d)', '_', name)

def tsv_to_mysql_and_export(tsv_file: str, db_config: dict, table_name: str, export_file: str):
    """
    将 TSV 文件内容导入 MySQL 数据库中的指定表格，并导出为压缩文件。

    :param tsv_file: TSV 文件路径
    :param db_config: 数据库连接配置字典 { 'host': 'localhost', 'user': 'root', 'password': 'password', 'database': 'dbname' }
    :param table_name: 目标表格的名称
    :param export_file: 导出的 SQL 文件路径
    """
    try:
        # 连接到 MySQL 数据库
        connection = mysql.connector.connect(
            host=db_config['host'],
            user=db_config['user'],
            password=db_config['password'],
            database=db_config['database']
        )
        cursor = connection.cursor()

        # 打开并读取 TSV 文件
        with open(tsv_file, 'r', encoding='utf-8') as file:
            tsv_reader = csv.reader(file, delimiter='\t')

            # 获取列名，并清洗列名中的特殊字符
            columns = next(tsv_reader)
            cleaned_columns = [clean_column_name(col) for col in columns]
            
            # 创建表格的 SQL 语句
            column_definitions = ', '.join([f'{col} VARCHAR(255)' for col in cleaned_columns])
            create_table_sql = f"CREATE TABLE IF NOT EXISTS {table_name} ({column_definitions});"
            cursor.execute(create_table_sql)

            # 创建 INSERT 语句
            placeholders = ', '.join(['%s'] * len(cleaned_columns))  # 为每一列生成一个 %s 占位符
            sql_insert = f"INSERT INTO {table_name} ({', '.join(cleaned_columns)}) VALUES ({placeholders})"

            # 插入数据
            for row in tsv_reader:
                cursor.execute(sql_insert, row)

            # 提交事务
            connection.commit()
            print(f"数据成功导入到 {table_name} 表。")

        # 导出 SQL 文件
        export_sql_file(cursor, table_name, export_file)

        # 压缩 SQL 文件
        compress_sql_to_zip(export_file)

    except mysql.connector.Error as err:
        print(f"数据库错误: {err}")
    
    finally:
        if cursor:
            cursor.close()
        if connection:
            connection.close()

def export_sql_file(cursor, table_name: str, export_file: str):
    """ 将数据库表导出为 SQL 文件 """
    with open(export_file, 'w', encoding='utf-8') as file:
        # 获取表的结构
        cursor.execute(f"SHOW CREATE TABLE {table_name}")
        create_table_sql = cursor.fetchone()[1]
        file.write(create_table_sql + ';\n')

        # 获取表数据并插入到 SQL 文件
        cursor.execute(f"SELECT * FROM {table_name}")
        rows = cursor.fetchall()

        for row in rows:
            insert_sql = f"INSERT INTO {table_name} VALUES ({', '.join(['%s'] * len(row))});"
            file.write(insert_sql % tuple(row) + '\n')

def compress_sql_to_zip(export_file: str):
    """ 压缩 SQL 文件到 ZIP 文件 """
    zip_file = export_file.replace('.sql', '.zip')
    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipf.write(export_file, os.path.basename(export_file))
    os.remove(export_file)  # 删除原始的 SQL 文件

    print(f"SQL 文件已压缩成 {zip_file}")

# 示例数据库配置
db_config = {
    'host': 'localhost',      # 数据库主机
    'user': 'root',           # 用户名
    'password': 'password',   # 密码
    'database': 'your_database'  # 数据库名称
}

# 使用函数进行转换
tsv_to_mysql_and_export('/mnt/ntc_data/wayne/Repositories/CRISPR/dual_20/NC_000024.10.tsv', db_config, 'your_table_name', 'exported_data.sql')
