def read_file(filename):
    import string
    f1=open(filename,'r')
    parts= []
    for line in f1:
        part = []
        for (word,i) in zip(string.split(line),range(8)):
            if i>=6:
                part.append(int(word))
            else:
                part.append(float(word))
        parts.append(part)
    return parts

def read_all_files(s):
    parts = []
    ok = True
    i = 100
    while ok:
        i+=1
        try:
            filename = s+'Parts00'+'%d'%i
            print filename
            parts += read_file(filename)
        except:
            ok = False
    return parts
    

old_parts = read_all_files('OLD')
new_parts = read_all_files('NEW')


ok = True
good_ones = 0
bad_ones = 0
for old_part in old_parts:
    idp = old_part[6]
    new_part = None
    for part in new_parts:
        if part[6]==idp:
            new_part = part
    if not old_part==new_part:
        ok = False
        print old_part
        print new_part
        bad_ones += 1
    else:
        good_ones += 1
if ok:
    print 'Success', good_ones
else:
    print 'Fail', good_ones, bad_ones
    

