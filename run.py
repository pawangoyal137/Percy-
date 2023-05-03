import subprocess
import threading
import re
import os
import signal

B = 32768
R = 32768

def extract_num(s):
    """
    Return the numbers in s
    """
    pattern = r"[-+]?\d*\.\d+|\d+"  # Regular expression for matching integers and floats
    match = re.search(pattern, s)
    if match:
        return float(match.group())
    else:
        return None

def read_server_output(id, server_process):
    total_time_ms = 0
    for line in server_process.stdout:
        if "Time to generate response" in line:
            total_time_ms+=extract_num(line)
        # print(f"Server{id}: {line}", end='') # process line here
    print(f"Server{id}: {total_time_ms}")
    
def start_server(id):
    print(f"Starting Server")
    from subprocess import Popen, PIPE, CalledProcessError
    server_process = subprocess.Popen(f"./pirserver big_database {id} {R} {B}", 
                                shell=True,
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT,
                                bufsize=1,
                                universal_newlines=True,
                                preexec_fn=os.setsid)
    
    for line in server_process.stdout:
        print(f"Server{id}: {line}", end='') # process line here
        if "Listening on port" in line:
            port_num = int(extract_num(line))
            print(f"Got port num = {port_num}")
            output_thread = threading.Thread(target=read_server_output, args=(id, server_process))
            break
    output_thread.start()
    return server_process, output_thread, port_num


def start_client(port_num1, port_num2):
    print(f"Starting Client")
    client_process = subprocess.Popen(f"./pirclient {R} {B} '1:localhost:{port_num1} 2:localhost:{port_num2}' 1 1", 
                                shell=True,
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT,
                                bufsize=1)
                                # universal_newlines=True)
    for line in client_process.stdout:
        print(f"Client: {line} \n", end='') # process line here

NUM_TRIALS = 2
server_process1, output_thread1, port_num1 = start_server(1)
server_process2, output_thread2, port_num2 = start_server(2)
for _ in range(NUM_TRIALS):
    start_client(port_num1, port_num2)
os.killpg(os.getpgid(server_process1.pid), signal.SIGTERM)
os.killpg(os.getpgid(server_process2.pid), signal.SIGTERM) 
output_thread1.join()
output_thread2.join()