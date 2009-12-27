import socket

PORT = 1233
CHUNK = 32768 * 8

def main():
    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    client_socket.connect(('localhost', PORT))

    for i in xrange(100000):
        client_socket.send('ASDFFFG ADFFFFG')
        a = client_socket.recv(CHUNK)

    client_socket.close()

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect(('localhost', PORT))
    print "DONE"

    for i in xrange(100):
        s.send("--matrix PAM250 --gap_open -10 ASDFFFG ADFFFFGZ")
        a = s.recv(CHUNK)
        #print a

    s.send('ASDF')
    print s.recv(CHUNK)
    s.send('ASDFFFG ADFFFFGZ')
    print s.recv(CHUNK)

    s.send('a')
    print s.recv(CHUNK)

    s.send('EXIT')

if __name__ == "__main__":
    main()
