import os
import argparse

def generate_full_directory_tree(current_path, output_file_handle, indent_prefix='', is_last_item=False):
    """
    递归生成目录树并写入文件句柄，修复了缩进和连接符问题。
    :param current_path: 当前正在处理的路径 (目录或文件)
    :param output_file_handle: 文件写入句柄
    :param indent_prefix: 当前层级的缩进前缀 (例如 '│   ' 或 '    ')
    :param is_last_item: 标识当前项是否是其父目录中的最后一个
    """
    # 获取当前目录下的所有文件和目录，并排序
    # 忽略隐藏文件和目录（以'.'开头），确保只处理可见文件
    items = [item for item in os.listdir(current_path) if not item.startswith('.')]
    items.sort()

    for i, item in enumerate(items):
        path = os.path.join(current_path, item)
        
        # 判断当前项是否是其父目录中的最后一个子项
        item_is_last = (i == len(items) - 1)
        
        # 根据是否是最后一个元素选择不同的连接符
        pointer = '└── ' if item_is_last else '├── '
        
        # 构建下一层级的缩进前缀
        # 如果当前项是最后一个，则后续子项的垂直线不再延伸
        next_indent = indent_prefix + ('    ' if item_is_last else '│   ')

        # 写入当前项
        if os.path.isdir(path):
            output_file_handle.write(f"{indent_prefix}{pointer}{item}/\n")
            # 递归调用处理子目录
            generate_full_directory_tree(path, output_file_handle, next_indent)
        else:
            output_file_handle.write(f"{indent_prefix}{pointer}{item}\n")

def merge_files_recursively(source_directory, output_file):
    """
    递归地将指定目录下及其所有子目录中文件的内容合并到一个文件中。
    输出文件将包含完整的目录结构列表，然后是每个文件的内容。
    """
    # 定义更长的分隔符
    delimiter = "################################################################################"

    try:
        with open(output_file, 'w', encoding='utf-8') as outfile:
            # --- 写入完整的目录结构和文件列表 ---
            outfile.write(f"{delimiter}\n")
            outfile.write(f"### FULL DIRECTORY AND FILE LIST START:\n")
            outfile.write(f"{delimiter}\n\n")

            # 获取根目录名
            root_dir_name = os.path.basename(os.path.normpath(source_directory))
            if not root_dir_name: # 如果是根目录，basename可能为空
                root_dir_name = source_directory

            # 根目录的显示
            # 特殊处理根目录，不带前置的连接符
            outfile.write(f"{root_dir_name}/\n")
            
            # 从根目录开始递归生成树，使用空字符串作为初始缩进
            generate_full_directory_tree(source_directory, outfile, '')

            outfile.write(f"\n{delimiter}\n")
            outfile.write(f"### FULL DIRECTORY AND FILE LIST END\n")
            outfile.write(f"{delimiter}\n\n\n")

            # --- 写入合并后的文件内容 ---
            outfile.write(f"{delimiter}\n")
            outfile.write(f"### MERGED FILE CONTENTS START:\n")
            outfile.write(f"{delimiter}\n\n")

            # 收集所有文件的路径，以便按排序顺序合并
            all_files_to_merge = []
            for root, _, files in os.walk(source_directory):
                # 过滤隐藏文件，确保只合并在目录树中显示的文件
                visible_files = [f for f in files if not f.startswith('.')]
                for filename in visible_files:
                    filepath = os.path.join(root, filename)
                    all_files_to_merge.append(filepath)
            
            # 对所有文件路径进行排序，以确保一致的合并顺序
            all_files_to_merge.sort()

            for filepath in all_files_to_merge:
                # 构建相对路径，使其在输出中更清晰
                relative_path = os.path.relpath(filepath, source_directory)

                outfile.write(f"{delimiter}\n")
                outfile.write(f"### FILE START: {relative_path}\n") # 文件开始标记，显示相对路径
                outfile.write(f"{delimiter}\n\n")
                try:
                    with open(filepath, 'r', encoding='utf-8') as infile:
                        outfile.write(infile.read())
                    outfile.write(f"\n\n{delimiter}\n")
                    outfile.write(f"### FILE END: {relative_path}\n") # 文件结束标记
                    outfile.write(f"{delimiter}\n\n") # 确保文件之间有额外的空行分隔
                except Exception as e:
                    print(f"Error reading file {relative_path}: {e}")
                    outfile.write(f"\n\n{delimiter}\n")
                    outfile.write(f"### ERROR READING FILE: {relative_path} - {e}\n")
                    outfile.write(f"{delimiter}\n\n")
            
            outfile.write(f"{delimiter}\n")
            outfile.write(f"### MERGED FILE CONTENTS END:\n")
            outfile.write(f"{delimiter}\n")

        print(f"所有文件已成功合并到：{output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

# --- 如何使用 ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="递归合并指定目录及其所有子目录下的文件。先列出完整目录结构，再合并文件内容，并使用长'#'符号作为标记。")
    parser.add_argument('source_dir', type=str,
                        help="包含要合并的文件的源目录路径。")
    parser.add_argument('output_file', type=str,
                        help="合并后的文件输出路径和名称。")

    args = parser.parse_args()

    # 检查源目录是否存在
    if not os.path.isdir(args.source_dir):
        print(f"错误：源目录 '{args.source_dir}' 不存在。请检查路径是否正确。")
    else:
        merge_files_recursively(args.source_dir, args.output_file)
