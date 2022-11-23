import numpy as np
import sys
import argparse
import pyslow5 as slow5
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.figsize'] = [12.0, 16.0]


'''
dRNA DNA barcode extraction.
looks for drop in signal and get's it as a segment

output structure:
fast5, readID, start, stop


This is to prototype read-time adapter detection.

Two methods:

1. Same as before, but we guess/use previous reads to build boundary conditions
   and compare each chunk against those until we hit the polyA

1a. We use variance to set the threshold for the end of the adapter detection
    as well as previous reads to set some initial boundary conditions to pick
    adapter state vs stall state at the start

2. Use variance to detect the polyA. Get mean of polyA region, then set boundary
   for adapter and then do adapter detection all at once. This has an advantage,
   where it might be faster and can just skip adapter detection and pull signal
   after polyA for readUntil

'''
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def get_chunks(signal, size=1000):
    """
    get chunks of signal as generator
    """
    batch = []
    for i in signal:
        batch.append(i)
        if len(batch) >= size:
            yield batch
            batch = []
    if len(batch) > 0:
        yield batch



def main():
    '''
    main function
    '''

    parser = MyParser(
        description="dRNA_segmenter - cut out adapter region of dRNA signal")
    #group = parser.add_mutually_exclusive_group()
    parser.add_argument("slow5",
                        help="slow5 file")
    # parser.add_argument("-w", "--window", type=int, default=250,
    #                     help="Window size")
    parser.add_argument("--chunksize", type=int, default=1200,
                        help="readUntil chunksize")
    parser.add_argument("-i", "--start_chunks", type=int, default=6,
                        help="num of chunks to store before processing")
    parser.add_argument("-c", "--corrector", type=int, default=1200,
                        help="Window size for increasing total error correction - better long segment detection")
    parser.add_argument("-w", "--window", type=int, default=300,    
                        help="Minimum segment window size to be detected")
    parser.add_argument("-e", "--error", type=int, default=5,
                        help="Allowable error in segment algorithm")
    parser.add_argument("-d", "--seg_dist", type=int, default=1800,
                        help="Maximum distance between 2 segments to be merged into 1")
    parser.add_argument("-m", "--min_seg_len", type=int, default=4000,
                        help="Minimum length of a segment to be constructed")
    parser.add_argument("-t", "--std_scale", type=float, default=0.9,
                        help="Scale factor of STDev about median")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="plot output")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    # arguments...put this into something better for Tansel
    # sig = args.signal         # signal file
    # SS = args.seq_sum          # Seq_sum.txt
    w = args.window      # window size
    # val_reads = sys.argv[3]
    # error = int(sys.argv[3])
    # thresh = int(sys.argv[4])
    s5 = slow5.Open(args.slow5, 'r')
    chunksize = args.chunksize
    start_chunks = [i for i in range(args.start_chunks)]
    reads = s5.seq_reads_multi(threads=8, batchsize=4000, pA=True, aux='all')
    for read in reads:

        readID = read['read_id']
        # sig = scale_outliers(read['signal'], args)
        sig = read['signal']
        # now feed algorithm with chunks of signal simulating real-time
        chunks = get_chunks(sig, size=chunksize)
        sig_store = np.array([])
        cum_sig = np.array([])
        sig_length = 0

            # this is the algo. Simple yet effective
        prev = False  # previous string
        err = 0       # total error
        prev_err = 0  # consecutive error
        c = 0         # counter
        w = args.corrector        # window to increase total error thresh
        seg_dist = args.seg_dist  # distance between 2 segs to be merged as one
        start = 0     # start pos
        end = 0       # end pos
        segs = []     # segments [(start, stop)]
        std_list = []
        med_list = []
        adapter_found = False
        for chunk_count, chunk in enumerate(chunks):
            # print("processing chunk: {}".format(chunk_count))
            # print("chunk {}".format(chunk))
            if chunk_count in start_chunks:
                # chunk_meds.append(median)
                # chunk_stdevs.append(stdev)
                sig_store = np.append(sig_store, chunk)
                median = np.median(chunk)
                stdev = np.std(chunk)
                std_list.append(stdev)
                med_list.append(median)
                # cum_sig = np.append(cum_sig, chunk)
                continue
            # parameter tuning visualisation
            # TODO: Put tuning plots here
            if len(sig_store) > 0:
                chunk = np.append(sig_store, chunk)
                sig_store = []
                # sig_length += len(chunk)
                # mn = sig.min()
                # mx = sig.max()
                mean = np.mean(chunk)
                # mean = 120
                median = np.median(chunk)
                # use this with outlier rejection to fix stdev thresholds
                stdev = np.std(chunk)
                # get variance
                var = np.var(chunk)
                top = median + (stdev * args.std_scale)
                # bot = median - (stdev * args.std_scale)
                # top = mean + (stdev * args.std_scale)
                # bot = mean - (stdev * args.std_scale)
            
            cum_sig = np.append(cum_sig, chunk)
           
            for i in range(len(chunk)):
                sig_length += 1
                a = chunk[i]
                if a < top: # If datapoint is within range
                    if not prev:
                        start = sig_length
                        prev = True
                    c += 1 # increase counter
                    w += 1 # increase window corrector count
                    if prev_err:
                        prev_err = 0
                    if c >= args.window and c >= w and not c % w: # if current window longer than detect limit, and corrector, and is divisible by corrector
                        err -= 1 # drop current error count by 1
                else:
                    if prev and err < args.error:
                        c += 1
                        err += 1
                        prev_err += 1
                        if c >= args.window and c >= w and not c % w:
                            err -= 1
                    elif prev and c >= args.window:
                        end = sig_length - prev_err # go back to where error stretch began for accurate cutting
                        prev = False
                        if segs and start - segs[-1][1] < seg_dist: # if segs very close, merge them
                            segs[-1][1] = end
                        else:
                            segs.append([start, end]) # save segment
                        c = 0
                        err = 0
                        prev_err = 0
                    elif prev:
                        prev = False
                        c = 0
                        err = 0
                        prev_err = 0
                    elif segs and (segs[-1][1]-segs[-1][0] >= args.min_seg_len) and sig_length - segs[-1][1] > seg_dist:
                        break_point = sig_length
                        # print("Break point: {}".format(break_point))
                        prev = False
                        c = 0
                        err = 0
                        prev_err = 0
                        adapter_found = True
                        break
                    else:
                        continue
            if adapter_found:
                break

        for a, b in segs:
            x, y = a, b
            print("{}\t{}\t{}".format(readID, x, y))
            break
            

    
        if args.plot:
            # print("plotting...")
            fig = plt.figure(1)
            #fig.subplots_adjust(hspace=0.1, wspace=0.01)
            ax = fig.add_subplot(111)
            # length = len(sig)
            # Show segment lines
            # readid, r_start, r_end =
            # i = int(r_start)
            # j = int(r_end)
            # if i > j:
            #     i, j = j, i
            # if i == j:
            #     j = j + 1
            # half_diff = (j-i) / 2
            # mid = i + half_diff
            # print(i, j, length)
            # print(sig[i:j])
            # y = max(sig[i:j]) + 20
            # ax.text(mid, y, readid)
            # ax.axvline(x=x, color='m')
            for i in [j for j in range(args.start_chunks+1)]:
                ax.axvline(x=i*1200, color='red')
            
            # show where it stopped processing more signal and broke out
            ax.axvline(x=break_point, color='blue')
            
            # print(std_list)
            for i in range(len(std_list)):
                st = std_list[i]
                med = med_list[i]
                xstart = i*1200
                xstop= (i+1)*1200
                ytop = med+st
                ybot = med-st
                # plt.axvspan(i*1200, (i+1)*1200, median-st, median+st, alpha=0.5, color='pink')
                # ax.axvspan(xstart, xstop, ymin=ybot, ymax=ytop, alpha=0.5, color='pink')
                # plt.axvspan(xstart, xstop, alpha=0.5, color='pink')
                ax.fill_between([xstart, xstop], ybot, ytop, alpha=0.8, color='violet')
            
            # ax.axvline(x=sig_length, color='blue')
            ax.axvspan(x, y, alpha=0.5, color='orange')
            ax.axhline(y=median, color='g')
            # ax.axhline(y=std, color='')
            ax.axhline(y=top, color='purple')
            plt.plot(sig, color='k')
            plt.show()
            # plt.savefig("test_{}.png".format(c))
            plt.clf()
            # sys.exit()




def scale_outliers(squig):
    ''' Scale outliers to within m stdevs of median '''
    k = (squig > 0) & (squig < 1200)
    return squig[k]


if __name__ == '__main__':
    main()
