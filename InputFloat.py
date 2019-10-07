#!/usr/bin/python
#
# this is a subroutine for userinput of a float number 
# inlcuding a default choice and timeout
# 
import sys
import time
import sys
if sys.platform=="win32":import msvcrt
else: import select, termios


def InputFloat( caption, default, timeout = 10):
    def writecaption(caption, default, input):
        #clear line from start
        sys.stdout.write('\r'+' '*79)
        #rewrite line from start
        sys.stdout.write('\r%s [%s]: %s'%(caption, default, input))
        sys.stdout.flush()
    
    try: default=float(default) # just check
    except: default =0.0        # and safeguard 
    
    if sys.platform=="win32":
        while msvcrt.kbhit(): msvcrt.getch()        #empty buffer  
        start_time = time.time()
        got_kbhit = False
        while True: #float coneversion loop
            input = ''    
            writecaption(caption, default, input)        
            while True: #read character from keyboard loop
                if msvcrt.kbhit():
                    got_kbhit = True
                    byte_arr = msvcrt.getche()
                    if ord(byte_arr) == 13: # enter_key
                        break
                    if ord(byte_arr) == 8:  # backspace
                        input = input[:-1]
                        writecaption(caption, default, input)
                    elif ord(byte_arr)>=32 and ord(byte_arr)<=126: #space .. "~" 
                        input += byte_arr.decode("unicode_escape") 
                    else: # clear trash caracters
                        while msvcrt.kbhit(): msvcrt.getch() #empty buffer
                        writecaption(caption, default, input)
                if not got_kbhit and (time.time() - start_time) > timeout:        
                    break
                time.sleep(0.01) # get CPU usage down
            if input=='': input=default
            try: 
                result=float(input); 
                break
            except:
                sys.stdout.write('\r'+' '*79)
                sys.stdout.write('\a\rInput Error: not a number')
                sys.stdout.flush()
                time.sleep(0.8)                                
        sys.stdout.write('\n'); sys.stdout.flush()  #move to next line
        while msvcrt.kbhit(): msvcrt.getch()        #empty buffer 
        return result        
        
    else: #linux
        termios.tcflush(sys.stdin, termios.TCIFLUSH) #empty buffer
        sys.stdout.write('\n\033[A') #scroll up one line        
        while True: #float coneversion loop
            input = '' 
            writecaption(caption, default, input)            
            i, o, e = select.select( [sys.stdin], [], [], timeout )
            if (i): 
               input = sys.stdin.readline().strip()
               sys.stdout.write('\033[A') #move cursor up
            if input=='': input=default; 
            try: 
                result=float(input); 
                break
            except:
                sys.stdout.write('\r'+' '*79)
                sys.stdout.write('\a\rInput Error: not a number')
                sys.stdout.flush()
                time.sleep(0.8)
        sys.stdout.write('\n'); sys.stdout.flush()   #move to next line
        termios.tcflush(sys.stdin, termios.TCIFLUSH) #empty buffer   
        return result 

#example usage
# 
#from InputFloat import InputFloat       
#answer = InputFloat('Input a number', 128.0, 30) 
#print( 'Result is %f' % answer) 








