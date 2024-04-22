import numpy as np
import sys
import argparse
import pyslow5 as slow5
import matplotlib.pyplot as plt
from matplotlib import rcParams
# rcParams['figure.figsize'] = [12.0, 16.0]
rcParams['figure.figsize'] = [20.0, 12.0]


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
    parser.add_argument("--chunksize", type=int, default=500,
                        help="readUntil chunksize")
    parser.add_argument("-i", "--start_chunks", type=int, default=4,
                        help="num of chunks to store before processing")
    parser.add_argument("-k", "--calc_chunk", type=int, default=2,
                        help="which chunk to start calculating threshold")
    parser.add_argument("-c", "--corrector", type=int, default=1200,
                        help="Window size for increasing total error correction - better long segment detection")
    parser.add_argument("-w", "--window", type=int, default=100,
                        help="Minimum segment window size to be detected")
    parser.add_argument("-e", "--error", type=int, default=5,
                        help="Allowable error in segment algorithm")
    parser.add_argument("-d", "--seg_dist", type=int, default=800,
                        help="Maximum distance between 2 segments to be merged into 1")
    parser.add_argument("-m", "--min_seg_len", type=int, default=3000,
                        help="Minimum length of a segment to be constructed")
    parser.add_argument("-n", "--no_err_thresh", type=int, default=800,
                        help="number of measurements to ignore errors")
    parser.add_argument("-t", "--std_scale", type=float, default=1.0,
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
    # w = args.window      # window size
    # val_reads = sys.argv[3]
    # error = int(sys.argv[3])
    # thresh = int(sys.argv[4])
    s5 = slow5.Open(args.slow5, 'r')
    chunksize = args.chunksize
    start_chunks = [i for i in range(args.start_chunks)]
    reads = s5.seq_reads_multi(threads=8, batchsize=4000, pA=True, aux='all')
    for read in reads:

        readID = read['read_id']
        # if readID != "0036c370-e1be-4b56-bc18-84eb9f4d4d41":
        #     continue
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
            if chunk_count in start_chunks:
                sig_store = np.append(sig_store, chunk)
                median = np.median(chunk)
                stdev = np.std(chunk)
                std_list.append(stdev)
                med_list.append(median)
                continue
            # parameter tuning visualisation
            # TODO: Put tuning plots here
            if len(sig_store) > 0:
                chunk = np.append(sig_store, chunk)
                sig_store = []
                mean = np.mean(chunk[args.calc_chunk:])
                # mean = 120
                median = np.median(chunk[args.calc_chunk:])
                # use this with outlier rejection to fix stdev thresholds
                stdev = np.std(chunk[args.calc_chunk:])
                # get variance
                var = np.var(chunk[args.calc_chunk:])
                top = median + (stdev * args.std_scale)

            
            cum_sig = np.append(cum_sig, chunk)
           
            for i in range(len(chunk)):
                sig_length += 1
                a = chunk[i]
                if a < top: # If datapoint is within range
                    if not prev:
                        start = sig_length
                        prev = True
                        err = 0
                    c += 1 # increase counter
                    if prev_err:
                        prev_err = 0
                    
                    # if current window longer than detect limit, and corrector, and is divisible by corrector
                    if c >= args.window and c >= w and not c % w and err > 0:
                        err -= 1 # drop current error count by 1

                else: # not within range
                    if prev and err < args.error: # winthin segment and less than error
                        c += 1
                        if sig_length >= args.no_err_thresh:
                            err += 1
                            prev_err += 1
                        if c >= args.window and c >= w and not c % w and err > 0:
                            err -= 1
                    elif prev: # within segment, above error, and greater than window
                        if c >= args.window:
                            end = sig_length - prev_err # go back to where error stretch began for accurate cutting
                            prev = False
                            if segs and start - segs[-1][1] < seg_dist: # if segs very close, merge them
                                segs[-1][1] = end
                            else:
                                segs.append([start, end]) # save segment
                            c = 0
                            err = 0
                            prev_err = 0
                        else: # within segment but not long enough
                            prev = False
                            c = 0
                            err = 0
                            prev_err = 0
                    elif segs and sig_length - segs[-1][1] > seg_dist:
                        break_point = sig_length
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
        
        if not adapter_found:
            break_point = sig_length

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
                ax.axvline(x=i*args.chunksize, color='red')
            
            # show where it stopped processing more signal and broke out
            ax.axvline(x=break_point, color='blue')
            
            # print(std_list)
            for i in range(len(std_list)):
                st = std_list[i]
                med = med_list[i]
                xstart = i*args.chunksize
                xstop= (i+1)*args.chunksize
                ytop = med+st
                ybot = med-st
                # plt.axvspan(i*1200, (i+1)*1200, median-st, median+st, alpha=0.5, color='pink')
                # ax.axvspan(xstart, xstop, ymin=ybot, ymax=ytop, alpha=0.5, color='pink')
                # plt.axvspan(xstart, xstop, alpha=0.5, color='pink')
                ax.fill_between([xstart, xstop], ybot, ytop, alpha=0.8, color='violet')
            
            # ax.axvline(x=sig_length, color='blue')
            ax.axvspan(x, y, alpha=0.5, color='orange')
            ax.axhline(y=median, color='g')
            # ax.axhline(y=90, color='r')
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
