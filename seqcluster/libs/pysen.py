
def pysenMMean (x, y):
        ymf = []
        total = 0
        for i in range(0, len(y)-1):
                total += y[i]
                ymf.append(total/(i+1))
        mf = 0
        xf = 0
        for i in range(0, len(y)-1):
                if ymf[i] > mf:
                        mf = ymf[i]
                        xf = i
        ymb = []
        total = 0
        for i in range(0, len(y)-1):
                ii = len(y)-1-i
                total += y[ii]
                ymb.append(total/(i+1))
        mb = 0
        xb = 0
        for i in range(0, len(y)-1):
                if ymb[i] > mb:
                        mb = ymb[i]
                        xb = i
        xb = len(y)-1 - xb
        return [x[xb], x[xf]]

# print "mmean %s" % pysenMMean(x , y)
