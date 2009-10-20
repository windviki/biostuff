import socket
from nwalign import global_align
import sys

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 1233
CHUNK = 32768
HOST = 'localhost'

server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.bind((HOST, PORT))
server.listen(5)
print "\nstarted server on %s:%i\n" % (HOST, PORT)


while True:
    client, address = server.accept()
    data = True 
    while data:
        try: 
            data = client.recv(CHUNK).rstrip()
            print "DATA:", data
            if data == "DONE":
                client.close()
                server.close()
                print "EXITING service"
                sys.exit(0)

            a, b = data.split(" ")
            r = global_align(a, b)
            client.send(" ".join(r))
        except Exception, e:
            try:
                client.send("ERROR:" + str(e))
            except socket.error:
                # they already closed...
                client.close()
                break


    client.close()
