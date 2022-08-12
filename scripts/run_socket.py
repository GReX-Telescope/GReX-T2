#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:20:32 2020

@author: gechen
"""

import socket
from astropy.io import ascii
import numpy as np
import pandas as pd

HOST = "127.0.0.1"
PORT = 12345

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((HOST, PORT))
s.listen(5)

while True:
    conn, add = s.accept()
    print(f"Connection from {add} has been established")
    # with conn:
    # ascii_letter = conn.recv(111111)
    ascii_letter = conn.recv(
        1
    )  # recieves an alphabet whose ASCII value is the size of the message

    if len(ascii_letter) > 0:

        # size = ord(ascii_letter.decode('utf-8'))      # ord() returns the ASCII value of a character
        # print(size)
        # print(ascii_letter)
        # print(ascii_letter.decode('utf-8'))

        # data = conn.recv(size)           # recieving the actual msg

        data = conn.recv(int(2e10))  # how to make sure that the size if larger enough?
        # data = data.split()

        # print(data)
        # print(type(data))
        print(len(data))
        data = conn.recv(len(data))

        candsfile = data.decode("utf-8")
        print(candsfile)
        # print(type(candsfile))
        # print(len(candsfile))

        # d = np.loadtxt(data)
        # print(d)
        # print(d.transpose())

        # tab = ascii.read(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'], guess=False, fast_reader=False, delimiter="\t")
        # print(tab)
        # print(tab['snr'])
