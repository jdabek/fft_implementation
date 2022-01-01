import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys

def readData(fin):
    trans = []
    with open(fin) as f:
        for line in f:
            l = line.strip()
            x = float(l.split(' ')[0])
            y = float(l.split(' ')[1])
            trans.append(x + y*1j)

    return np.array(trans)

def writeData(fout, data):
    f = open(fout, 'w')
    for d in data:
        f.write(format(d.real, '.20f') + ' ' + format(d.imag, '.20f'))
        f.write('\n')
    f.close()

def callFftImplementation(specs):
    subprocess.run(['./fft_implementation.o', specs['trans'], specs['fin'], specs['fout']])
    dataout = readData(specs['fout'])
    return dataout

def performValidation(specs, datain=None):
    t = np.arange(specs['nsamp'])

    if datain == None:
        datain = np.random.randn(specs['nsamp']) + np.random.randn(specs['nsamp'])*1j
    writeData(specs['fin'], datain)

    if specs['trans'] == 'fft':
        trans1 = np.fft.fft(datain)
    else:
        trans1 = np.fft.ifft(datain)
    trans2 = callFftImplementation(specs)

    rmserr = np.sqrt(np.average(np.abs(trans2 - trans1)**2))
    # Normalize the fft to be like ifft:
    if specs['trans'] == 'fft':
        rmserr /= specs['nsamp']
    specs['rmserr'] = rmserr

    plt.clf()
    plt.plot(t, trans1.real, 'c-', t, trans1.imag, 'y-', linewidth=5.0)
    plt.plot(t, trans2.real, 'b--', t, trans2.imag, 'm--')
    plt.legend(['xpython', 'ypython', 'ximpl', 'yimpl'])
    plt.savefig(specs['fimg'])

    return datain

def processValidations(specs, sampCounts):
    allSpecs = []
    accuracyFft = []
    accuracyIfft = []
    for i in range(len(sampCounts)):
        specs['nsamp'] = sampCounts[i]
        minmax = '(from: ' + str(sampCounts[0]) + ' to: ' + str(sampCounts[-1]) + ')'
        print('Processing validation nsamp ' + minmax + ': ' + str(sampCounts[i]))
        specs['trans'] = 'fft'
        datain = performValidation(specs)
        allSpecs.append(specs)
        accuracyFft.append([specs['nsamp'], specs['rmserr']])
        specs['trans'] = 'ifft'
        datain = performValidation(specs)
        allSpecs.append(specs)
        accuracyIfft.append([specs['nsamp'], specs['rmserr']])

    return [np.array(accuracyFft), np.array(accuracyIfft), allSpecs]

if __name__ == '__main__':
    args = sys.argv

    if len(args) != 3:
        print('Please give 2 integer arguments to validate the fft implmementation.')
        print('If given e.g. "1 130", the fft and ifft are validated')
        print('from 1 sample to 130 samples with increments of 1,')
        print('with unit-normally distributed random noise.')
        sys.exit('Incompatible arguments')

    nsampMin = int(args[1])
    nsampMax = int(args[2])
    if nsampMin < 1:
        nsampMin = 1
    if nsampMin > nsampMax:
        sys.exit('Incompatible arguments: ' + str(args[1:]))

    specs = {}

    specs['fin'] = 'input.txt'
    specs['fout'] = 'output.txt'
    specs['fimg'] = 'output.png'

    sampCounts = [i for i in range(nsampMin, nsampMax+1)]

    [accuracyFft, accuracyIfft, allSpecs] = processValidations(specs, sampCounts)

    plt.clf()
    plt.plot(accuracyFft[:,0], accuracyFft[:,1], 'b-')
    plt.plot(accuracyIfft[:,0], accuracyIfft[:,1], 'r-')
    plt.legend(['fft', 'ifft'])
    plt.title('Accuracies with unit-normally distributed noise\n(fft normalized: divided by nsamp)')
    plt.xlabel('nsamp')
    plt.ylabel('rmserr')
    plt.savefig('rmserrs.png')

    isFailed = False
    threshold = 10e-15
    for i in range(len(accuracyFft)):
        if abs(accuracyFft[i][1]) > threshold:
            isFailed = True
            outstr = 'Accuracy fft FAILED at nsamp '
            outstr += str(int(accuracyFft[i][0]))
            outstr += ' w/rmserr '
            outstr += str(accuracyFft[i][1])
            print(outstr)
        if abs(accuracyIfft[i][1]) > threshold:
            isFailed = True
            outstr = 'Accuracy ifft FAILED at nsamp '
            outstr += str(int(accuracyIfft[i][0]))
            outstr += ' w/rmserr '
            outstr += str(accuracyIfft[i][1])
            print(outstr)
    if isFailed:
        print('Validation tests FAILED w/rmserr threshold ' + str(threshold))
    else:
        print('Validation tests PASSED w/rmserr threshold ' + str(threshold))
