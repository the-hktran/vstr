def AdjustLine(line, N):
    p = len(line) - 1
    adjusted_line = ''
    for i in range(p):
        adjusted_line += str(int(line[i]) - N) + ' '
    adjusted_line += line[p]
    adjusted_line = str(p) + ' ' + adjusted_line + '\n'
    return adjusted_line

def CFOURConverter(Outname, N):
    o = open(Outname, 'w')
    f = open("cubic", 'r')

    line = f.readline().split()
    newline = AdjustLine(line, N)
    o.write(newline)

    while True:
        line = f.readline().split()
        if not line:
            break;
        newline = AdjustLine(line, N)
        o.write(newline)

    f = open("quartic", 'r')
    while True:
        line = f.readline().split()
        if not line:
            break;
        newline = AdjustLine(line, N)
        o.write(newline)

if __name__ == '__main__':
    CFOURConverter('jf_input.inp', 7)

