import re
import argparse
import os

keywords = ["input", "output", "log", "jobid","reason","threads","resources","["]
def is_line_not_start_with_keywords_and_not_empty(line, keywords):
    # 判断line是否为空行或以keywords中任意一个字符串开头
    return not (line.strip() == "" or any(line.strip().startswith(keyword) for keyword in keywords))

# 测试数据，这里以列表形式提供了特定的字符串


def extract_command_lines(log_file):
    rule_commands = {}
    current_rule = None
    with open(log_file, 'r') as file:
        rule_start = False
        
        for line in file:
            #print(rule_commands)
            line = line.strip()
            # Check for rule header lines
            match = re.match(r'^rule (\w+):', line)
            if rule_start == False:
                
                if match:
                    current_rule = match.group(1)
                    rule_commands[current_rule] = []
                    rule_start = True
                    
            elif rule_start == True:
                #print(line)
                if match:
                    current_rule = match.group(1)
                    rule_commands[current_rule] = []
            #Check for command lines under rules
                else:
                    if is_line_not_start_with_keywords_and_not_empty(line,keywords):
                        rule_commands[current_rule].append(line)
                #print(rule_commands)
    return rule_commands


def save_command_lines(rule_commands,cmdfile,outpath):
    cmdallfile = open(cmdfile,"a")
    for rule, commands in rule_commands.items():
        cmdallfile.write("\n".join(commands))
        
        cmdrulefile = os.path.join(outpath ,rule+".sh")
        with open(cmdrulefile, 'w') as file:
            file.write('\n\n'.join(commands))
    cmdallfile.close()


if __name__ == "__main__":
    #print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())) + "annotation of fusion ")
    parser = argparse.ArgumentParser(description='get_cmd')
    parser.add_argument('-f','--logfile', help="log file with commands")
    parser.add_argument('-o','--outpath',help="outpath of results")
    parser.add_argument('-l','--cmdfile',help="cmdfile")
    
    args = parser.parse_args()
    log_file = args.logfile
    outpath = args.outpath
    cmdfile = args.cmdfile
    
    print(log_file)
    rule_commands = extract_command_lines(log_file)
    #print("check",rule_commands)
    save_command_lines(rule_commands,cmdfile,os.path.join(outpath,"cmd"))