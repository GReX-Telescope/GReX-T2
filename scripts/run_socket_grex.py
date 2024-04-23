import sys

import socket
from grex_t2 import socket_grex
import logging as logger

logger.basicConfig(filename="output.log", encoding="utf-8", level=logger.DEBUG)

HOST = "127.0.0.1"
PORT = 12345

def main(trigger=True, args=sys.argv):
    min_dm = float(args[1])
    max_ibox = float(args[2])
    min_snr = float(args[3])
    min_snr_t2out = float(args[4])
    
    # Use roughly 8 seconds as a gulp size
    gulpsize = 16384 * 8

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    # Create a UDP socket
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    # Ensure that you can reconnect
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    # Bind the socket to the port
    server_address = (HOST, PORT)
    s.bind(server_address)

    print("Connected to socket %s:%d. Triggering set to %s" % (HOST, PORT, trigger))
    logger.info("Connected to socket %s:%d. Triggering set to %s" % (HOST, PORT, trigger))

    candsfile = ["", "", "", "", ""]

    last_trigger_time = 0.0
    # Outer loop that runs as long as T2 is running
    while True:
        candstr_list = ''
        cand_count = 0
        # Inner loop for chunks of Heimdall output data
        while True:
            # Recieve 512 bytes
            data, address = s.recvfrom(512)

            # Waiting for end of text. When chunk is done, break inner loop
            if len(data)==1 and data==b'\x03':
                break
            
            candstr = data.decode("utf-8")
            
            # Removing the \n from the end of line
            candstr_list += candstr
            cand_count += 1

        print("Number of candidates %d" % cand_count)

        if cand_count > 0:
            print("Filtering", "last trig was ", last_trigger_time)
            last_trigger_time = socket_grex.filter_candidates(candstr_list, trigger=trigger, last_trigger_time=last_trigger_time, min_dm=min_dm, max_ibox=max_ibox, min_snr=min_snr, min_snr_t2out=min_snr_t2out)
            print("Finished filtering")
            
        continue

