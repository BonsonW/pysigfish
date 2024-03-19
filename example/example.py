#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pysigfish
import pyslow5



class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


class RUread:
    '''
    batch class to mimic readUntil object
    '''

    def __init__(self, id, number, chunk_length, raw_data):
        # print(id, number, chunk_length, raw_data)
        self.id = id
        self.number = number
        self.chunk_length = chunk_length
        self.raw_data = np.ndarray.tobytes(raw_data)


def real_main2(args):
    '''
    test
    '''
    channels = args.channels

    s5 = pyslow5.Open(args.slow5, 'r')
    recs = s5.seq_reads_multi(aux='all', threads=args.threads, batchsize=channels, pA=True)
    records = []
    for rec in recs:
        records.append(rec)

    # pysig = pysigfish.start(args.reference, args.paf, channels=channels, threads=args.threads, DEBUG=1)
    pysig = pysigfish.start(args.reference, args.paf, channels=channels, threads=1, DEBUG=1)
    batch = []
    dddtype = 'f'
    round = 0
    F = open('decisions.tsv', 'w')
    print("round: {}".format(round))
    C = 0
    G = 0
    while round < 10:
        for rec in records:
            read = RUread(str(G), G, rec["len_raw_signal"], rec["signal"])
            channel = C + 1
            batch.append([channel, read])
            C += 1
            G += 1
            # print("channel: {}".format(channel))
            if C == 512:
                status = pysig.process_batch(batch, dddtype)
                for ch in status.keys():
                    read_number = status[ch][1]
                    read_id = status[ch][2]
                    decision = status[ch][3]
                    print(read_id, ch, read_number, decision, sep='\t', file=F)
                C = 0
                round += 1
                batch = []
                print("round: {}".format(round))
                break
    print("done!")
    F.close()
    s5.close()


def real_main3(args):
    '''
    test
    '''
    CHANNELS = 10
    CHUNK_SIZE = 1200
    ROUNDS=50

    s5 = pyslow5.Open(args.slow5, 'r')
    recs = s5.seq_reads_multi(aux='all', threads=args.threads, batchsize=CHANNELS, pA=True)
    records = []
    for rec in recs:
        records.append(rec)

    pysig = pysigfish.start(args.reference, args.paf, channels=CHANNELS, threads=args.threads, DEBUG=1)
    batch = []
    dddtype = 'f'
    F = open('decisions.tsv', 'w')
    for round in range(ROUNDS):
        print("round: {}".format(round))
        for channel in range(CHANNELS):
            # print(records[channel+1]["read_id"])
            print("chunks in read {} {} chunks:{} total:{}".format(records[channel+1]["read_id"], channel+1, (round+1)*CHUNK_SIZE, records[channel+1]['len_raw_signal']))
            if (round+1)*CHUNK_SIZE >= records[channel+1]['len_raw_signal']:
                print("No more chunks in read {} {} chunks:{} total:{}".format(records[channel+1]["read_id"], channel+1, (round+1)*CHUNK_SIZE, records[channel+1]['len_raw_signal']))
                continue
            read = RUread(records[channel+1]["read_id"], 0, CHUNK_SIZE, records[channel+1]["signal"][CHUNK_SIZE*round:(CHUNK_SIZE*(round+1))+1])
            batch.append([channel+1, read])
            # print("channel: {}".format(channel))
        # for i, j in batch:
            # print(j.id, i, j.number)
        status = pysig.process_batch(batch, dddtype)

        for ch in status.keys():
            read_number = status[ch][1]
            read_id = status[ch][2]
            decision = status[ch][3]
            # print(read_id, ch, read_number, decision, sep='\t', file=F)
            print(read_id, ch, read_number, decision, sep='\t')
        batch = []
    print("done!")
    F.close()
    s5.close()



def main():

    parser = MyParser(description="test pysigfish",
    epilog="Citation:...",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-c", "--channels", type=int, default=512,
                        help="number of channels")
    parser.add_argument("-s", "--slow5", default="/home/jamfer/Data/SK/hasindu/sequin_reads.blow5",
                        help="slow5 file")
    parser.add_argument("-r", "--reference", default="/home/jamfer/Data/cheat_MinKNOW/gencode.v40.10p.ref.fa",
                        help="reference fasta")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="num of threads to use")
    parser.add_argument("-p", "--paf", default="./out.paf",
                        help="paf output")
    parser.add_argument("-d", "--dev",type=int, default=1,
                        help="log level, 0 is off, 1 is on")


    args = parser.parse_args()

    # if len(sys.argv) == 1:
    #     parser.print_help(sys.stderr)
    #     sys.exit(1)

    # real_main2(args)
    real_main3(args)


if __name__ == '__main__':
    main()
    

    '''
    int real_main2(const char *slow5file, const char *fasta_file, int num_thread, int8_t full_ref){

    slow5_file_t *sp = slow5_open(slow5file,"r");
    if(sp==NULL){
       fprintf(stderr,"Error in opening file\n");
       exit(EXIT_FAILURE);
    }
    slow5_rec_t *rec = NULL;
    int ret=0;

    int channels=512;

    sigfish_opt_t opt;
    opt.num_thread = num_thread;
    opt.debug_paf = "-";
    opt.no_full_ref = !full_ref;
    sigfish_state_t *state = init_sigfish(fasta_file, channels, opt);
    sigfish_read_t reads[channels];

    int round=0;
    while (ret==0){
        int i = 0;
        while (i < channels && (ret = slow5_get_next(&rec,sp))>= 0) {
            reads[i].channel = i+1;
            reads[i].read_number = round;
            reads[i].len_raw_signal = rec->len_raw_signal;
            reads[i].read_id = (char *)malloc(strlen(rec->read_id)+1);
            strcpy(reads[i].read_id,rec->read_id);
            assert(reads[i].read_id != NULL);
            reads[i].raw_signal = (float*)malloc(sizeof(float)*rec->len_raw_signal);
            for (int j=0; j<rec->len_raw_signal; j++){
                reads[i].raw_signal[j] = TO_PICOAMPS(rec->raw_signal[j],rec->digitisation,rec->offset,rec->range);
            }
            i++;
        }
        fprintf(stderr,"round %d: %d reads loaded\n",round,i);

        enum sigfish_status *status = process_sigfish(state, reads, i);
        for(int j=0;j<i;j++){
            fprintf(stderr,"channel %d: %d\n",j+1,status[j]);
        }
        free(status);
        for(int j=0; j<i; j++){
            free(reads[j].raw_signal);
            free(reads[j].read_id);
        }
        fprintf(stderr,"round %d: %d reads processed\n",round,i);
        round++;

    }

    free_sigfish(state);
    slow5_rec_free(rec);
    slow5_close(sp);

    return 0;

    }
    '''
