import socket

PORT = 1233
CHUNK = 32768 * 8

"""
client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
client_socket.connect(('localhost', PORT))

for i in xrange(10000):
    client_socket.send('ASDFFFG ADFFFFG')
    a = client_socket.recv(CHUNK)

client_socket.close()
"""

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect(('localhost', PORT))
print "DONE"

for i in xrange(10000):
    s.send('--matrix PAM250 --gap 10 ASDFFFG ADFFFFGZ')
    a = s.recv(CHUNK)
    print a


s.send('ASDF')
print s.recv(CHUNK)
s.send('ASDFFFG ADFFFFGZ')
print s.recv(CHUNK)

s.send('a')
print s.recv(CHUNK)

s.send('EXIT')
