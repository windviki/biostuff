import socket
from nwalign import global_align

PORT = 1233
CHUNK = 32768
HOST = 'localhost'

server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.bind((HOST, PORT))
server.listen(5)
print "started server on %s:%i" % (HOST, PORT)


while True:
    client, address = server.accept()
    data = True 
    while data:
        try: 
            data = client.recv(CHUNK).strip()
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
