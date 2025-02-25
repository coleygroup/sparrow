import subprocess

def run_commands(cmds):
    for cmd in cmds: 
        print(cmd)
        subprocess.call(cmd, shell=True)

def load_commands(file: str = 'examples/automoldesigner/scripts/cmd_list.txt'): 
    with open(file, 'r') as f:
        cmds = f.read().splitlines()
    return cmds 

if __name__=='__main__': 
    cmds = load_commands()
    run_commands(cmds)
