import dnaio
from xopen import xopen
import os
import traceback
import io
from multiprocessing import Process, Pipe, Queue
from typing import BinaryIO
from qualtrim_new import quality_trim_index
from ultraplex.modifiers import AdapterCutter, ModificationInfo
from ultraplex.adapters import BackAdapter
import gzip
import glob
import time
import argparse
import shutil
from pathlib import Path
import logging
import sys
from math import log10, floor
import numpy as np

# @profile
def hamming_distance(seq1, seq2, limit=0):
    """
    Returns the Hamming distance between equal-length sequences.
    :param seq1: first sequence. - str
    :param seq2: second sequence. - str
    :param limit: max distance limit before aborting (returning limit + 1).
    :return: the edit distance.
    """
    i = 0
    sum = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            sum += 1
        if limit > 0 and sum > limit:
            return limit + 1
    return sum

def subglobal_distance(s2, s1):
    """
    Computes the edit distance for a sub-global alignment
    of a sequence s2 against a sequence s1.
    Mismatches and indels both score as 1. Overhanging parts of s1 do not count.
    :param s1: the longer (probe) sequence.
    :param s2: the shorter sought sequence.
    :return: the minimum edit distance
    """
    xLen = len(s1)
    yLen = len(s2)

    if xLen < yLen:
        raise ValueError("Sub-global edit distance is undefined for sequences "
                         "where the probe is shorter than the aligned sequence.")

    d = np.empty([xLen + 1, yLen + 1], dtype=np.uint32)

    # Initialize array
    d[:, 0] = 0
    d[0, :] = np.arange(yLen + 1)

    # Perform DP.
    for x in range(1, xLen + 1):
        # Fill matrix.
        for y in range(1, yLen + 1):
            d[x, y] = min(d[x - 1, y] + 1, d[x, y - 1] + 1, d[x - 1, y - 1] + int(s1[x - 1] != s2[y - 1]))

    # Find min for sub-global alignment so that all of s2 is covered,
    # but not necessarily all of s1 sequence.
    mini = 1000000
    iPos = 0
    i = xLen
    while i > 0:
        if d[i, yLen] < mini:
            mini = d[i, yLen]
            iPos = i
        i -= 1

    # Return min distance
    return mini

def user_trim(read, trim_sequences):
    """Simply helper function to remove
    helper sequences from a barcode"""
    
    if len(read.sequence) < trim_sequences[-1][-1]:
        raise ValueError("Invalid trimming sequences given " \
                         "The number of positions given must be even and they must fit into the barcode length.")

    prev_start = 0
    prev_end = 0
    for start,end in trim_sequences:
        offset = prev_end - prev_start
        read.sequence = read.sequence[:(start-offset)] + read.sequence[(end-offset):]
        read.qualities = read.qualities[:(start-offset)] + read.qualities[(end-offset):]
        prev_start = start
        prev_end = end
    return read

def round_sig(x, sig=2):
    try:
        to_return = round(x, sig - int(floor(log10(abs(x)))) - 1)
    except:
        to_return = 0
    return to_return


def make_5p_bc_dict(barcodes, min_score, dont_build_reference):
    """
	this function generates a dictionary that matches each possible sequence
	from the read with the best 5' barcode
	"""
    max_edit_distance = 8
    if dont_build_reference:
        return {"dont_build": True}
    else:
        first_bc = barcodes[0]
        seq_length = len(first_bc.replace("N", ""))

        # check all the barcodes are the same length (ignoring UMIs)
        # for bc in barcodes:
        #     assert len(bc.replace("N", "")) == seq_length, "Your experimental barcodes are different lengths."

        seqs = make_all_seqs(seq_length)
        
        # print(seqs)
        # print("")

        # trim sequences to desired length
        # create dictionary
        barcode_dictionary = {}

        for seq in seqs:
            barcode_dictionary[seq] = score_barcode_for_dict(seq, barcodes, max_edit_distance)
        return barcode_dictionary

# @profile
def remove_Ns_from_barcodes(barcodes):
    # returns a dictionary with keys of the N-removed barcode, and values of the original barcode
    # barcodes_no_N = {a.replace("N", ""): a for a in barcodes}
    # barcodes_no_N["no_match"] = "no_match"
    return barcodes

# # @profile
# def score_barcode_for_dict(seq, barcodes, max_edit_distance, Ns_removed=False):
#     """
# 	this function scores a given sequence against all the barcodes. It returns the winner with Ns included.
# 	"""

#     max_edit_distance=8
#     if not Ns_removed:
#         barcodes = remove_Ns_from_barcodes(barcodes)
#     barcodes_no_N = barcodes.keys()

#     if seq in barcodes_no_N:  # no need to check all barcodes
#         winner = barcodes[seq]  # barcode WITH Ns included
#     # elif min_score == len(barcodes_no_N):  # i.e. no matches allowed, and seq not in barcodes
#     #     winner = "no_match"
#         # print(seq)
#         # print(winner)
#         # print("")
#     else:  # mismatches allowed so need to check
#         dists = {}

#         for this_bc in barcodes_no_N:
#             # # score the barcode against the read, penalty for N in the read
#             # score = sum(a == b for a, b in zip(this_bc, seq))
#             # scores[this_bc] = score
#             if this_bc != "no_match":
#                 #dist = subglobal_distance(this_bc, seq)
#                 dist = hamming_distance(this_bc, seq)
#                 dists[this_bc] = dist

#         # Find the best score
#         min_dist = min(dists.values())

#         if min_dist > max_edit_distance:
#             winner = "no_match"
#         else:
#             # check that there is only one barcode with the max score
#             filtered = [a for a, b in dists.items() if b == min_dist]
#             # if len(filtered) > 1:
#             #     winner = "no_match-ambiguous" #need to decide how to fix the multiple matches situation
#             # else:  # if there is only one
#             #     winner = barcodes[filtered[0]]  # barcode WITH Ns included
#             winner = barcodes[filtered[0]]  # barcode WITH Ns included
#     #     print(min_dist)
#     # print(seq)
#     # print(winner)
#     # print("")
#     return winner

# @profile

def score_barcode_for_dict(seq, barcodes, max_edit_distance, Ns_removed):
    """
	this function scores a given sequence against all the barcodes. It returns the winner with Ns included.
	"""

    if seq in barcodes:  # no need to check all barcodes
        winner = seq  # barcode WITH Ns included
    # elif min_score == len(barcodes_no_N):  # i.e. no matches allowed, and seq not in barcodes
    #     winner = "no_match"
        # print(seq)
        # print(winner)
        # print("")
    else:  # mismatches allowed so need to check
        dists = {}

        for this_bc in barcodes:
            # # score the barcode against the read, penalty for N in the read
            # score = sum(a == b for a, b in zip(this_bc, seq))
            # scores[this_bc] = score
            if this_bc != "no_match":
                #dist = subglobal_distance(this_bc, seq)
                dist = hamming_distance(this_bc, seq)
                dists[this_bc] = dist

        # Find the best score
        min_dist = min(dists.values())

        if min_dist > max_edit_distance:
            winner = "no_match"
        else:
            # check that there is only one barcode with the max score
            filtered = [a for a, b in dists.items() if b == min_dist]
            # if len(filtered) > 1:
            #     winner = "no_match-ambiguous" #need to decide how to fix the multiple matches situation
            # else:  # if there is only one
            #     winner = barcodes[filtered[0]]  # barcode WITH Ns included
            #print(dists)
            winner = filtered[0]  # barcode WITH Ns included
    #     print(min_dist)
    # print(seq)
    # print(winner)
    # print("")
    return winner

class ReaderProcess(Process):
    """
	Read chunks of FASTA or FASTQ data (single-end or paired) and send to a worker.

	The reader repeatedly

	- reads a chunk from the file(s)
	- reads a worker index from the Queue
	- sends the chunk to connections[index]

	and finally sends the stop token -1 ("poison pills") to all connections.
	"""

    # @profile
    def __init__(self, file, connections, queue, buffer_size):#, i2):
        # /# Setup the reader process

        """
		queue -- a Queue of worker indices. A worker writes its own index into this
			queue to notify the reader that it is ready to receive more data.
		connections -- a list of Connection objects, one for each worker.
		"""
        super().__init__()
        self.file = file
        self.connections = connections
        self.queue = queue
        self.buffer_size = buffer_size
        #self.file2 = i2

    # @profile
    def run(self):
        # if self.stdin_fd != -1:
        # 	sys.stdin.close()
        # 	sys.stdin = os.fdopen(self.stdin_fd) #/# not sure why I need this!

        try:
            with xopen(self.file, 'rb') as f:
                for chunk_index, chunk in enumerate(dnaio.read_chunks(f, self.buffer_size)):
                    #print("Ahsley needs help",chunk_index,sys.getsizeof(chunk))
                    self.send_to_worker(chunk_index, chunk)

            # Send poison pills to all workers
            for _ in range(len(self.connections)):
                worker_index = self.queue.get()
                self.connections[worker_index].send(-1)
        except Exception as e:
            # TODO better send this to a common "something went wrong" Queue
            for connection in self.connections:
                connection.send(-2)
                connection.send((e, traceback.format_exc()))

    #@profile
    def send_to_worker(self, chunk_index, chunk, chunk2=None):
        worker_index = self.queue.get()  # get a worker that needs work
        connection = self.connections[worker_index]  # find the connection to this worker
        connection.send(chunk_index)  # send the index of this chunk to the worker
        connection.send_bytes(chunk)  # /# send the actual data to this worker
        if chunk2 is not None:
            connection.send_bytes(chunk2)


class InputFiles:
    """
	this is from cutadapt - basically just creates a dnaio object
	"""

    def __init__(self, file1: BinaryIO, interleaved: bool = False):
        self.file1 = file1
        self.interleaved = interleaved

    def open(self):
        return dnaio.open(self.file1,
                          interleaved=self.interleaved, mode="r")


def make_all_seqs(l):
    """
	Makes all possible sequences, including Ns, of length l
	"""
    if l > 8:
        print("Warning - large barcodes detected!")
        print("It may be faster to use option '--dont_build_reference'!")

    nts = ['A', "C", "G", "T", "N"]

    all_seqs = nts

    for i in range(l - 1):
        new_seqs = []
        for seq in all_seqs:
            for nt in nts:
                new_seqs.append(seq + nt)
        all_seqs = new_seqs

    return (all_seqs)


def rev_c(seq):
    """
    simple function that reverse complements a given sequence
    """
    tab = str.maketrans("ACTGN", "TGACN")
    # first reverse the sequence
    seq = seq[::-1]
    # and then complement
    seq = seq.translate(tab)
    return seq


class WorkerProcess(Process):  # /# have to have "Process" here to enable worker.start()
    """
	The worker repeatedly reads chunks of data from the read_pipe, runs the pipeline on it
	and sends the processed chunks to the write_pipe.

	To notify the reader process that it wants data, processes then writes out.
	"""

    # @profile
    def __init__(self, index,
                 read_pipe, need_work_queue,
                 output_directory,
                 barcodes, coordinates,
                 save_name,
                 total_demultiplexed,
                 total_reads_assigned,
                 ultra_mode,
                 min_score_5_p, 
                 ignore_no_match,
                 dont_build_reference,
                 trim_sequences,
                 kmers_dict):
        
        super().__init__()
        self._id = index  # the worker id
        self._read_pipe = read_pipe  # the pipe the reader reads data from
        self._need_work_queue = need_work_queue  # worker adds its id to this queue when it needs work

        self._total_demultiplexed = total_demultiplexed  # a queue which keeps track of the total number of reads processed
        self._total_reads_assigned = total_reads_assigned  # a queue which keeps track of the total number of reads assigned to sample files

        self._save_name = save_name  # the name to save the output fastqs
        #self._five_p_barcodes_pos, self._five_p_umi_poses = find_bc_and_umi_pos(barcodes)
        self._five_p_bc_dict = make_5p_bc_dict(barcodes, min_score_5_p, dont_build_reference)
        self._min_score_5_p = min_score_5_p  
        self._ultra_mode = False
        self._output_directory = output_directory

        self._ignore_no_match = ignore_no_match
        self._barcodes_no_N = remove_Ns_from_barcodes(barcodes)

        self._trim_sequences = trim_sequences
        self._kmers_dict = kmers_dict
        self._coordinates = coordinates

    # @profile
    def run(self):
        
        
        while True:  # /# once spawned, this keeps running forever, until poison pill recieved
            # Notify reader that we need data
            self._need_work_queue.put(self._id)

            # /# get some data
            chunk_index = self._read_pipe.recv()

            # /# check there's no error
            if chunk_index == -1:  # /# poison pill from Sina
                # reader is done
                break
            elif chunk_index == -2:
                # An exception has occurred in the reader
                e, tb_str = self._read_pipe.recv()
                raise e

            # /# otherwise if we have no error, run...
            # /# get some bytes
            data = self._read_pipe.recv_bytes()
            infiles = io.BytesIO(data)

            # /# process the reads
            #processed_reads = []
            this_buffer_list = []
            reads_written = 0
            assigned_reads = 0


            for read in InputFiles(infiles).open():
                reads_written += 1
                #umi = ""
                
                read = user_trim(read, [(10,24),(34,50)])
                # trim_sequences = [(30, 34), (80, 88), (90, 100)]
                trim_sequences = self._trim_sequences
                if trim_sequences:
                    read = user_trim(read, trim_sequences)
                
                
                # /# demultiplex at the 5' end ###
                read.name = read.name.replace(" ", "").replace("/", "").replace("\\",
                                                                                "")  # remove bad characters

                #print(read)
                read, five_p_bc = five_p_demulti(read,
                                                                    #self._five_p_barcodes_pos,
                                                                    #self._five_p_umi_poses,
                                                                    #self._five_p_bc_dict,
                                                                    #add_umi=True,
                                                                    barcodes_no_N=self._barcodes_no_N,
                                                                    min_score=self._min_score_5_p,
                                                                    kmers_dict=self._kmers_dict)
                                                                    #keep_barcode=self._keep_barcode)
                                                                    
                
                
                

                if five_p_bc != "no_match":
                    x_coord = self._coordinates[five_p_bc][0]
                    y_coord = self._coordinates[five_p_bc][1]
                    read.name = read.name + " B0:Z:" + five_p_bc + " B1:Z:" + x_coord + " B2:Z:" + y_coord
                    this_buffer_list.append(read)
                    #print(read.name)
                    
                    # filename = 'matched-4000.fastq.gz'

                    # if os.path.exists(filename):
                    #     append_write = 'ab'  # append if already exists
                    # else:
                    #     append_write = 'wb'  # make a new file if not

                    # with gzip.open(filename, append_write) as file:
                    #     this_out = []
                    #     this_out.append("@" + read.name)
                    #     this_out.append(read.sequence)
                    #     this_out.append("+")
                    #     this_out.append(read.qualities)
                    #     output = '\n'.join(this_out) + '\n'
                    #     file.write(output.encode())

            # ## Write out! ##
            # for read in this_buffer_list:
            #     write_temp_files(output_dir=self._output_directory,
            #                     save_name=self._save_name,
            #                     # demulti_type=demulti_type,
            #                     worker_id=self._id,
            #                     read=read)
            #                     # ultra_mode=self._ultra_mode,
            #                     # ignore_no_match=self._ignore_no_match)
            ## Write out! ##
            # for read in this_buffer_list:
            #     write_temp_files(output_dir=self._output_directory,
            #                     save_name=self._save_name,
            #                     #demulti_type=demulti_type,
            #                     worker_id=self._id,
            #                     read=read)
            #                     #ultra_mode=self._ultra_mode,
            #                     #ignore_no_match=self._ignore_no_match)
            

            write_temp_files(output_dir=self._output_directory,
                            save_name=self._save_name,
                            #demulti_type=demulti_type,
                            worker_id=self._id,
                            reads=this_buffer_list)
                            #ultra_mode=self._ultra_mode,
                            #ignore_no_match=self._ignore_no_match)
            
            # LOG reads processed
            prev_total = self._total_demultiplexed.get()
            new_total = prev_total[0] + reads_written
            if new_total - prev_total[1] >= 1_000_000:
                prog_msg = str(new_total // 1_000_000) + ' million reads processed'
                print(prog_msg)
                logging.info(prog_msg)
                last_printed = new_total
            else:
                last_printed = prev_total[1]
            self._total_demultiplexed.put([new_total, last_printed])

            # LOG reads assigned
            prev_total = self._total_reads_assigned.get()
            new_total = prev_total + assigned_reads
            self._total_reads_assigned.put(new_total)
# @profile
def write_temp_files(output_dir, save_name, worker_id, reads):
    # write_this = True  # assume true
    # if "no_match" in demulti_type and ignore_no_match:
    #     write_this = False

    # if ultra_mode and write_this:
    #     # /# work out this filename
    #     filename = output_dir + 'ultraplex_' + save_name + demulti_type + '_tmp_thread_' + str(
    #         worker_id) + '.fastq'

    #     if os.path.exists(filename):
    #         append_write = 'a'  # append if already exists
    #     else:
    #         append_write = 'w'  # make a new file if not

    #     with open(filename, append_write) as file:
    #         this_out = []
    #         for counter, read in enumerate(reads):
    #             # Quality control:
    #             assert len(read.name.split("rbc:")) <= 2, "Multiple UMIs in header!"

    #             if counter == 0:
    #                 umi_l = len(read.name.split("rbc:")[1])
    #             assert len(read.name.split("rbc:")[1]) == umi_l, "UMIs are different lengths"
    #             ## combine into a single list
    #             this_out.append("@" + read.name)
    #             this_out.append(read.sequence)
    #             this_out.append("+")
    #             this_out.append(read.qualities)

    #         output = '\n'.join(this_out) + '\n'
    #         # print(output)
    #         file.write(output)

    # elif write_this:
        # /# work out this filename
        # filename = output_dir + 'ultraplex_' + save_name + demulti_type + '_tmp_thread_' + str(
        #     worker_id) + '.fastq.gz'

    filename = output_dir + 'ultraplex_' + save_name + 'tmp_thread_' + str(
        worker_id) + '.fastq'

    if os.path.exists(filename):
        append_write = 'a'  # append if already exists
    else:
        append_write = 'wb'  # make a new file if not

    with gzip.open(filename, append_write) as file:
        this_out = []
        # for counter, read in enumerate(reads):

        #     # Quality control:
        #     assert len(read.name.split("rbc:")) <= 2, "Multiple UMIs in header!"

        #     if counter == 0:
        #         umi_l = len(read.name.split("rbc:")[1])
        #     assert len(read.name.split("rbc:")[1]) == umi_l, "UMIs are different lengths"
        #     ## combine into a single list
        for read in reads:
            this_out.append("@" + read.name)
            this_out.append(read.sequence)
            this_out.append("+")
            this_out.append(read.qualities)

        output = '\n'.join(this_out) + '\n'
        file.write(output.encode())
# @profile
def write_tmp_files(output_dir, save_name, demulti_type, worker_id, reads,
                    ultra_mode, ignore_no_match):
    write_this = True  # assume true
    if "no_match" in demulti_type and ignore_no_match:
        write_this = False

    if ultra_mode and write_this:
        # /# work out this filename
        filename = output_dir + 'ultraplex_' + save_name + demulti_type + '_tmp_thread_' + str(
            worker_id) + '.fastq'

        if os.path.exists(filename):
            append_write = 'a'  # append if already exists
        else:
            append_write = 'w'  # make a new file if not

        with open(filename, append_write) as file:
            this_out = []
            for counter, read in enumerate(reads):
                # Quality control:
                assert len(read.name.split("rbc:")) <= 2, "Multiple UMIs in header!"

                if counter == 0:
                    umi_l = len(read.name.split("rbc:")[1])
                assert len(read.name.split("rbc:")[1]) == umi_l, "UMIs are different lengths"
                ## combine into a single list
                this_out.append("@" + read.name)
                this_out.append(read.sequence)
                this_out.append("+")
                this_out.append(read.qualities)

            output = '\n'.join(this_out) + '\n'
            # print(output)
            file.write(output)

    elif write_this:
        # /# work out this filename
        filename = output_dir + 'ultraplex_' + save_name + demulti_type + '_tmp_thread_' + str(
            worker_id) + '.fastq.gz'

        if os.path.exists(filename):
            append_write = 'ab'  # append if already exists
        else:
            append_write = 'wb'  # make a new file if not

        with gzip.open(filename, append_write) as file:
            this_out = []
            for counter, read in enumerate(reads):

                # Quality control:
                assert len(read.name.split("rbc:")) <= 2, "Multiple UMIs in header!"

                if counter == 0:
                    umi_l = len(read.name.split("rbc:")[1])
                assert len(read.name.split("rbc:")[1]) == umi_l, "UMIs are different lengths"
                ## combine into a single list
                this_out.append("@" + read.name)
                this_out.append(read.sequence)
                this_out.append("+")
                this_out.append(read.qualities)

            output = '\n'.join(this_out) + '\n'
            # print(output)
            file.write(output.encode())

#@profile
def get_kmers_dict(barcodes,k):
    kmers_dict = {}
    
    for barcode in barcodes:
        prev = ''
        kmer = 0
        offset = 0

        width = len(barcode)
        
        for i in range(width):
            x = barcode[i]
            # print("i: "+str(i))
            # print("x: "+x)
            # print("prev: "+prev)
            # print("kmer: "+str(kmer))

            if i==0 or x==prev: #if starting or staying on same kmer
                if kmer < k:
                    kmer+=1

            else: #if ending a kmer
                if kmer>1:
                    key = (kmer, prev, offset)
                    #print(str(kmer) + " " + prev + " " + str(offset))
                    if key in kmers_dict: kmers_dict[key].extend([barcode])
                    else: kmers_dict[key] = [barcode]
                    offset = i
                kmer=1 
            if i == width-1:
                if kmer>1:
                    key = (kmer, x, offset)
                    #print(str(kmer) + " " + prev + " " + str(offset))
                    if key in kmers_dict: kmers_dict[key].extend([barcode])
                    else: kmers_dict[key] = [barcode]
                    offset = i
            prev=x
            #print("")

    #first_three_items = list(kmers_dict.items())[:3]
    #print(first_three_items)
    return kmers_dict
            

    # first_three_items = list(kmers_dict.items())[:3]
    # print(first_three_items)
    # return kmers_dict

#A_3[1] = [tgatgtctcccttttagcttttaaaa, tgcttaggggd]
#A_3[1] = [tgatgtctcccttttagcttttaaaa, tgcttaggggd]


    # A = [""]*k
    # G = [""]*k
    # T = [""]*k
    # C = [""]*k
    
    # for barcode in barcodes:
    #     prev = barcode[0]
    #     kmer = 0
    #     for x, i in barcode:
    #         if x==prev:
    #             if kmer < k:
    #                 kmer+1
    #         else:
    #             if x == "A": A[kmer] = {barcode:i}
    #             if x == "G": A[kmer] = {barcode:i}
    #             if x == "T": A[kmer] = {barcode:i}
    #             if x == "C": A[kmer] = {barcode:i}
                
    #             prev=x
    #             kmer=0
    
# def get_kmers(barcode,k):
#     kmers = []
#     prev = barcode[0]
#     kmer = 0
#     start = False
    
#     for x, i in barcode:
#         if x == 'N': 
#             start = True
#             prev = barcode[i+1]
    
#         if start == True:
#             if x==prev:
#                 if kmer < k:
#                     kmer+1
#             else:
#                 key = prev + "_" + str(kmer) + str(i)
#                 kmers.append(key)
#                 prev=x
#                 kmer=0
#     return kmers

#@profile
def get_kmers(sequence,k):
    kmers = []
    prev = ''
    kmer = 0
    offset = 0
    width = len(sequence)
        
    for i in range(width):
        x = sequence[i]

        if i==0 or x==prev: #if starting or staying on same kmer
            if kmer < k:
                kmer+=1
        else: #if ending a kmer
            if kmer > 1:
                key = (kmer, prev, offset)
                if prev !='N': kmers.append(key)
                offset = i
            kmer=1 
        if i == width-1: # if end of the string
            if kmer > 1:
                key = (kmer, x, offset)
                if x !='N': kmers.append(key)
                offset = i
        prev=x
    return kmers

    # for x, i in barcode:
    #     if prev == 'N':
    #         prev = x
    #     if x != 'N': 
    #         prev = barcode[i]
    #         if x==prev:
    #             if kmer < k:
    #                 kmer+1
    #         else:
    #             key = prev + "_" + str(kmer) + str(i)
    #             kmers.append(key)
    #             prev=x
    #             kmer=0
    # return kmers

    # for barcode in barcodes:
    #     prev = ''
    #     kmer = 0
    #     offset = 0

    #     width = len(barcode)
        
    #     for i in range(width):
    #         x = barcode[i]
    #         # print("i: "+str(i))
    #         # print("x: "+x)
    #         # print("prev: "+prev)
    #         # print("kmer: "+str(kmer))

    #         if i==0 or x==prev: #if starting or staying on same kmer
    #             if kmer < k:
    #                 kmer+=1

    #         else: #if ending a kmer
    #             key = (str(kmer)+prev, offset)
    #             # print("key:" + str(kmer)+prev+ ", " + str(offset))
    #             if key in kmers_dict: kmers_dict[key].extend([barcode])
    #             else: kmers_dict[key] = [barcode]
    #             offset = i
    #             kmer=1 
    #         if i == width-1:
    #             key = (str(kmer)+x, offset)
    #             #print("key:" + str(kmer)+prev+ ", " + str(offset))
    #             if key in kmers_dict: kmers_dict[key].extend([barcode])
    #             else: kmers_dict[key] = [barcode]
    #             offset = i
    #         prev=x
    #         #print("")

    # # first_three_items = list(kmers_dict.items())[:3]
    # # print(first_three_items)
    # return kmers_dict

# start = time.time()
#@profile
def get_candidates(sequence, kmers_dict):
    """
    Returns candidate barcodes for a read barcode
    as a list of barcodes
    First, split the barcode in kmers and then
    finds all the candidate barcodes of the kmers
    :param read_barcode the barcode from which to get candidates
    :return a list of candidates barcodes
    """
    # NOTE probably faster to keep kmer_offsets in memory as we will call
    #      this function several times with the same barcode but we get a penalty in memory use

    k=2
    sequence_kmers = get_kmers(sequence,k)
    pool = {}
    for key in sequence_kmers:
        kmer_len = key[0]
        nt = key[1]
        pos = key[2]
        
        if key in kmers_dict:
            pool.update({key: kmers_dict[key]})
        
        up_pos = pos
        down_pos = pos
        
        
        while up_pos < k + pos:
            new_key = (kmer_len,nt,up_pos)
            if new_key in kmers_dict and new_key not in pool:
                pool.update({new_key: kmers_dict[new_key]})
            new_key = (kmer_len+1,nt,up_pos)
            if new_key in kmers_dict and new_key not in pool:
                pool.update({new_key: kmers_dict[new_key]})
            up_pos+=1

        thresh = max(0, down_pos-k)
        while down_pos >= thresh:
            new_key = (kmer_len,nt,down_pos)
            if new_key in kmers_dict and new_key not in pool:
                pool.update({new_key: kmers_dict[new_key]})
            new_key = (kmer_len+1,nt,down_pos)
            if new_key in kmers_dict and new_key not in pool:
                pool.update({new_key: kmers_dict[new_key]})
            down_pos-=1
        
        unique_barcode_set = set()
        for lst in pool.values():
            unique_barcode_set.update(lst)
        # end = time.time()
        # print(end - start)
        #print(len(unique_barcode_set))
    return list(unique_barcode_set)

    
    # # Iterate all the kmer-offset combinations found in the input barcode
    # for kmer, offset in kmers_offsets:
    #     # Obtain all the barcodes that matched for the current kmer
    #     try:
    #         hits = kmer2seq[kmer]
    #     except KeyError:
    #         continue
    #     # For each true barcode containing read's kmer.
    #     # Hit refers to barcode and hit_offsets to where the kmer was in the barcode
    #     for hit, hit_offsets in list(hits.items()):
    #         # if no_offset_speedup:
    #         #     # NON-OPTIMIZED CASE
    #         #     # For each kmer in read (typically incremented by k positions at a time).
    #         #     candidates[hit] = 0
    #         #     continue
    #         # OPTIMIZED CASE
    #         # For each kmer in read (typically incremented by k positions at a time).
    #         min_penalty = 100000000
    #         # For each position that kmer occurred in the true barcode.
    #         # Get the minimum penalty
    #         for hit_offset in hit_offsets:
    #             # Kmer may be shifted overhang positions without penalty, due to subglobal alignment.
    #             penalty = max(0, abs(offset - hit_offset) - pre_overhang - post_overhang)
    #             if penalty < min_penalty: min_penalty = penalty
    #         # Assign the min penalty to the candidate (if exists already take max)
    #         # TODO if there are several equal barcode candidates for different kmers,
    #         #      why keep the max penalty and not an average?
    #         candidates[hit] = max(min_penalty, candidates[hit])
            
    # # Clear out all candidates with a forced offset penalty greater than the max edit distance and return
    # return [hit for hit,penal in list(candidates.items()) if penal <= max_edit_distance]

# def get_candidates(sequence):
#     """
#     Returns candidate barcodes for a read barcode
#     as a list of barcodes
#     First, split the barcode in kmers and then
#     finds all the candidate barcodes of the kmers
#     :param read_barcode the barcode from which to get candidates
#     :return a list of candidates barcodes
#     """
#     # NOTE probably faster to keep kmer_offsets in memory as we will call
#     #      this function several times with the same barcode but we get a penalty in memory use

#     kmers_offsets = ku.get_kmers(sequence, k)

    
#     # Iterate all the kmer-offset combinations found in the input barcode
#     for kmer, offset in kmers_offsets:
#         # Obtain all the barcodes that matched for the current kmer
#         try:
#             hits = kmer2seq[kmer]
#         except KeyError:
#             continue
#         # For each true barcode containing read's kmer.
#         # Hit refers to barcode and hit_offsets to where the kmer was in the barcode
#         for hit, hit_offsets in list(hits.items()):
#             # if no_offset_speedup:
#             #     # NON-OPTIMIZED CASE
#             #     # For each kmer in read (typically incremented by k positions at a time).
#             #     candidates[hit] = 0
#             #     continue
#             # OPTIMIZED CASE
#             # For each kmer in read (typically incremented by k positions at a time).
#             min_penalty = 100000000
#             # For each position that kmer occurred in the true barcode.
#             # Get the minimum penalty
#             for hit_offset in hit_offsets:
#                 # Kmer may be shifted overhang positions without penalty, due to subglobal alignment.
#                 penalty = max(0, abs(offset - hit_offset) - pre_overhang - post_overhang)
#                 if penalty < min_penalty: min_penalty = penalty
#             # Assign the min penalty to the candidate (if exists already take max)
#             # TODO if there are several equal barcode candidates for different kmers,
#             #      why keep the max penalty and not an average?
#             candidates[hit] = max(min_penalty, candidates[hit])
            
#     # Clear out all candidates with a forced offset penalty greater than the max edit distance and return
#     return [hit for hit,penal in list(candidates.items()) if penal <= max_edit_distance]

# @profile
def five_p_demulti(read, kmers_dict, barcodes_no_N=[], min_score=0):
    """
    this function demultiplexes on the 5' end
    """
    sequence_length = len(read.sequence)
    candidates = get_candidates(read.sequence, kmers_dict)
    #print(read.sequence)
    winner = score_barcode_for_dict(read.sequence, candidates, min_score, Ns_removed=True)
    if sequence_length < len(winner):  # read is too short to contain barcode
        winner = "no_match"
    # if sequence_length > max(five_p_bc_pos):
    #     # find best barcode match
    #     # this_bc_seq = ''.join([read.sequence[i] for i in five_p_bc_pos])
    #     this_bc_seq = read.sequence
    #     if "dont_build" in five_p_bc_dict:
    #         # print(this_bc_seq)
    #         winner = score_barcode_for_dict(this_bc_seq, barcodes_no_N, min_score, Ns_removed=True)
    #     else:
    #         winner = five_p_bc_dict[this_bc_seq]

    #     # store what sequence will be removed
    #     if sequence_length < len(winner):  # read is too short to contain barcode
    #         winner = "no_match"
            
    return read, winner

# @profile
def find_bc_and_umi_pos(barcodes):
    """
	This function finds the coordinates of the umi and barcodes nucleotides
	"""
    bcs_poses = {}
    umi_poses = {}
    for bc in barcodes:
        bcs_poses[bc] = [i for i in range(len(bc)) if bc[i] != "N"]
        umi_poses[bc] = [i for i in range(len(bc)) if bc[i] == "N"]

    # /# we assume that the barcode is always the same
    bc_pos = bcs_poses[barcodes[0]]

    umi_poses["no_match"] = []

    return bc_pos, umi_poses

# @profile
def start_workers(n_workers, input_file, need_work_queue, #adapter,
                  barcodes, coordinates, save_name, total_demultiplexed, total_reads_assigned,
                  min_score_5_p, #three_p_mismatches, linked_bcds, three_p_trim_q,
                  ultra_mode, output_directory, #final_min_length, q5, i2, adapter2, min_trim,
                  ignore_no_match, dont_build_reference, #keep_barcode, trim_sequences):
                  trim_sequences,
                  kmers_dict,buffer_size):
    """
	This function starts all the workers
	"""
    workers = []
    all_conn_r = []
    all_conn_w = []

    total_demultiplexed.put([0, 0])  # [total written, last time it was printed] - initialise [0,0]
    total_reads_assigned.put(0)

    for index in range(n_workers):
        # create a pipe to send data to this worker
        conn_r, conn_w = Pipe(duplex=False)
        all_conn_r.append(conn_r)
        all_conn_w.append(conn_w)

        worker = WorkerProcess(index=index,
                               read_pipe=conn_r,  # this is the "read_pipe"
                               need_work_queue=need_work_queue,  # worker tells the reader it needs work
                               output_directory=output_directory,
                               barcodes=barcodes,
                               coordinates=coordinates,
                               save_name=save_name,
                               total_demultiplexed=total_demultiplexed,
                               total_reads_assigned=total_reads_assigned,
                               ultra_mode=ultra_mode,
                               min_score_5_p=min_score_5_p,
                               ignore_no_match=ignore_no_match,
                               dont_build_reference=dont_build_reference,
                               trim_sequences=trim_sequences,
                               kmers_dict=kmers_dict
                               )

        worker.start()
        workers.append(worker)

    return workers, all_conn_r, all_conn_w

# @profile
def concatenate_files(save_name, ultra_mode,
                      sbatch_compression,
                      output_directory,
                      compression_threads=8):
    """
	this function concatenates all the files produced by the 
	different workers, then sends an sbatch command to compress
	them all to fastqs.
	"""
    # First, file all the unique file names we have, ignoring threads
    all_names = glob.glob(output_directory + "ultraplex_" + save_name + '*')

    all_types = []  # ignoring threads
    for name in all_names:
        this_type = name.split("_tmp_thread_")[0]

        if this_type not in all_types:
            all_types.append(this_type)  # this_type contains directory if applicable

    # now concatenate them
    if ultra_mode:
        for this_type in all_types:
            # find all files with this barcode (or barcode combination)
            filenames = sorted(glob.glob(this_type + '*'))  # this type already has output directory
            # then concatenate
            command = ''
            for name in filenames:
                command = command + name + ' '
            command = 'cat ' + command + ' > ' + this_type + '.fastq'
            os.system(command)

            for name in filenames:
                os.remove(name)

            print("Compressing with pigz...")
            c_thread_n = '-p' + str(compression_threads)
            if sbatch_compression:
                os.system(
                    'sbatch -J compression -c ' + str(
                        compression_threads) + ' --time 4:00:00 --wrap="pigz ' + c_thread_n + ' ' + this_type + '.fastq"')
            else:
                os.system('pigz ' + c_thread_n + ' ' + this_type + '.fastq')

        # check if compression is complete
        if sbatch_compression:
            finished = False
            print("Compressing....")
            while not finished:
                # assume it's complete
                complete = True
                # now actually check
                for this_type in all_types:
                    filename = glob.glob(this_type + '*')

                    if '.gz' not in filename[0]:
                        complete = False

                if complete:
                    finished = True
                    print("Compression complete!")
                else:
                    time.sleep(1)
    else:  # if not ultra_mode
        for this_type in all_types:
            # find all files with this barcode (or barcode combination)
            filenames = sorted(glob.glob(this_type + '*'))
            # then concatenate
            command = ''
            for name in filenames:
                command = command + name + ' '
            command = 'cat ' + command + ' > ' + this_type + '.fastq.gz'
            os.system(command)

            for name in filenames:
                os.remove(name)

# @profile
def clean_files(output_directory, save_name):
    files = glob.glob(output_directory + 'ultraplex_' + save_name + '*')
    for file in files:
        os.remove(file)

# @profile
def process_bcs(tsv, mismatch_5p):
    barcodes = []
    coordinates = {}

    counter = 0  # all 5' barcodes must be consistent

    with open(tsv, 'r') as file:
        for row in file:
            counter += 1
            # First, find if theres a comma
            line = row.split('\t')
            coordinates[line[0]]=(line[1],line[2].strip()) #in case \n
            
    barcodes =  list(coordinates.keys())
    match_5p = len(barcodes[0]) - mismatch_5p

    return barcodes, coordinates, match_5p

# @profile
def check_enough_space(output_directory, input_file,
                       ignore_space_warning, ultra_mode):
    # First, find the free space on the output directory
    if output_directory == "":
        output_directory = os.getcwd()
    total, used, free = shutil.disk_usage(output_directory)

    if ultra_mode:
        multiplier = 0.098
    else:
        multiplier = 0.98

    # Find the size of the input file
    input_file_size = Path(input_file).stat().st_size
    if ignore_space_warning:
        if not input_file_size < multiplier * free:
            print("WARNING! System may not have enough free space to demultiplex")
            print("(Warning has been ignored)")
    else:
        assert input_file_size < free * multiplier, "Not enough free space. To ignore this warning use option --ignore_space_warning"

# @profile
def check_N_position(bcds, type):
    # checks that UMI positions result in consistent barcode
    if not len(bcds) == 0:
        for counter, bcd in enumerate(bcds):
            # find positions of non-N
            non_N = [a for a, b in enumerate(bcd) if b != "N"]

            if type == "5":
                # then look for first non N
                ref_pos = min(non_N)
            else:
                # look for last non_n
                ref_pos = len(bcd) - max(
                    non_N)  # not just max(non_N) because need to allow for different UMI length at 5' end of 3' bcd

            if counter == 0:
                correct_pos = ref_pos
            else:
                assert ref_pos == correct_pos, "UMI positions not consistent"

# @profile
# def main(buffer_size=int(4 * 1024 ** 2)):  # 4 MB
def main(buffer_size=int(50 * 1024 ** 1)):  # 4 MB
    start = time.time()

    ## PARSE COMMAND LINE ARGUMENTS ##

    parser = argparse.ArgumentParser(description='Ultra-fast demultiplexing of fastq files.')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    # input
    required.add_argument('-i', "--inputfastq", type=str, required=True,
                          help='fastq file to be demultiplexed')
    # barcodes csv
    required.add_argument('-b', "--barcodes", type=str, required=True,
                          help='barcodes for demultiplexing in tsv format')
    # output directory
    optional.add_argument('-d', "--directory", type=str, default="", nargs='?',
                          help="optional output directory")
    optional.add_argument('-bs', "--buffersize", type=int, default=50, nargs='?',
                          help="Buffer Size in kb [DEFAULT 50]")
    # 5' mismatches
    optional.add_argument('-m5', "--fiveprimemismatches", type=int, default=1, nargs='?',
                          help='number of mismatches allowed for 5prime barcode [DEFAULT 1]')

    # threads
    optional.add_argument('-t', "--threads", type=int, default=4, nargs='?',
                          help='threads [DEFAULT 4]')

    # name of output file
    optional.add_argument('-o', "--outputprefix", type=str, default="demux", nargs='?',
                          help='prefix for output sequences [DEFAULT demux]')
    # use sbatch compression in ultra mode
    optional.add_argument('-sb', "--sbatchcompression", action='store_true', default=False,
                          help='whether to compress output fastq using SLURM sbatch')
    # ultra mode
    optional.add_argument('-u', "--ultra", action='store_true', default=False,
                          help='whether to use ultra mode, which is faster but makes very large temporary files')
    # free space ignore warning
    optional.add_argument('-ig', "--ignore_space_warning", action='store_true', default=False,
                          help='whether to ignore warnings that there is not enough free space')
    optional.add_argument("-inm", "--ignore_no_match", action="store_true", default=False,
                           help="Do not write reads for which there is no match.")
    optional.add_argument("-dbr", "--dont_build_reference", default=False, action="store_true",
                          help="Skip the reference building step - for long barcodes this will improve speed.")
    optional.add_argument("-ts","--trim-sequences", type=int, default=None, nargs='+', 
                        help="Trims from the barcodes in the input file\n" \
                        "The bases given in the list of tuples as START END START END .. where\n" \
                        "START is the integer position of the first base (0 based) and END is the integer\n" \
                        "position of the last base.\nTrimmng sequences can be given several times.")

    parser._action_groups.append(optional)
    args = parser.parse_args()

    output_directory = args.directory

    if not output_directory == "":
        if not output_directory[len(output_directory) - 1] == "/":
            output_directory = output_directory + "/"

    # Make output directory
    if not output_directory == "":
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

    logging.basicConfig(level=logging.DEBUG, filename=output_directory + "ultraplex_" + str(start) + ".log",
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")

    print(args)
    logging.info(args)

    file_name = args.inputfastq
    barcodes_tsv = args.barcodes
    mismatch_5p = args.fiveprimemismatches
    threads = args.threads
    save_name = args.outputprefix
    sbatch_compression = args.sbatchcompression
    ultra_mode = args.ultra
    ignore_space_warning = args.ignore_space_warning
    ignore_no_match = args.ignore_no_match
    dont_build_reference = args.dont_build_reference
    buffer_size = args.buffersize * 1024

    if ultra_mode:
        print("Warning - ultra mode selected. This will generate very large temporary files!")

    if not ultra_mode:
        if sbatch_compression:
            print("sbatch_compression can only be used in conjunction with ultra mode")
            print("setting sbatch_compression to false")
            sbatch_compression = False
    
    if args.trim_sequences is not None \
    and (len(args.trim_sequences) % 2 != 0 or min(args.trim_sequences) < 0): 
        raise ValueError("Invalid trimming sequences given " \
                         "The number of positions given must be even and they must fit into the barcode length.")
    
    # Make the input trim coordinates a list of tuples
    trim_sequences = None
    if args.trim_sequences is not None:
        trim_sequences = list()
        sorted_seqs = args.trim_sequences
        #sorted_seqs.sort() #removed for now because of -1 testing

        for i in range(len(sorted_seqs) - 1):
            if i % 2 == 0:
                trim_sequences.append((sorted_seqs[i], 
                                       sorted_seqs[i+1]))

    # assert output_directory=="" or output_directory[len(output_directory)-1]=="/", "Error! Directory must end with '/'"

    #check_enough_space(output_directory, file_name, ignore_space_warning, ultra_mode, i2)
    check_enough_space(output_directory, file_name, ignore_space_warning, ultra_mode)

    # process the barcodes csv
    #five_p_bcs, three_p_bcs, linked_bcds, min_score_5_p, sample_names = process_bcs(barcodes_tsv, mismatch_5p)
    barcodes, coordinates, min_score_5_p = process_bcs(barcodes_tsv, mismatch_5p)
    
    
    # start = time.time()
    kmers_dict = get_kmers_dict(barcodes, len(barcodes[0]))
    # end = time.time()
    # print(end - start)
    #print(kmers_dict)
    
    #check_N_position(barcodes, "5")  # check 3' later so that different 5' barcodes can have different types of 3' bcd

    # remove files from previous runs
    clean_files(output_directory, save_name)

    # /# Make a queue to which workers that need work will add
    # /# a signal
    need_work_queue = Queue()
    total_demultiplexed = Queue()
    total_reads_assigned = Queue()


    # /# make a bunch of workers
    workers, all_conn_r, all_conn_w = start_workers(n_workers=threads,
                                                    input_file=file_name,
                                                    need_work_queue=need_work_queue,
                                                    barcodes=barcodes,
                                                    coordinates=coordinates,
                                                    save_name=save_name,
                                                    total_demultiplexed=total_demultiplexed,
                                                    total_reads_assigned=total_reads_assigned,
                                                    min_score_5_p=min_score_5_p,
                                                    ultra_mode=ultra_mode,
                                                    output_directory=output_directory,
                                                    ignore_no_match=ignore_no_match,
                                                    dont_build_reference=dont_build_reference,
                                                    trim_sequences=trim_sequences,
                                                    kmers_dict=kmers_dict,
                                                    buffer_size=buffer_size)
    
    print("Demultiplexing...")
    reader_process = ReaderProcess(file_name, all_conn_w,
                                   need_work_queue, buffer_size)
    reader_process.daemon = True
    reader_process.run()

    # concatenate_files(save_name, ultra_mode, sbatch_compression, output_directory, sample_names)
    concatenate_files(save_name, ultra_mode, sbatch_compression, output_directory)
    total_processed_reads = total_demultiplexed.get()[0]
    runtime_seconds = str((time.time() - start) // 1)
    finishing_msg = "Demultiplexing complete! " + str(
        total_processed_reads) + ' reads processed in ' + runtime_seconds + ' seconds'
    print(finishing_msg)
    logging.info(finishing_msg)

    # More stats for logging
    total_ass = total_reads_assigned.get()
    total_ass_percent = (total_ass / total_processed_reads) * 100
    assmsg = str(total_ass) + " (" + str(
        round_sig(total_ass_percent, 3)) + "%) reads correctly assigned to sample files"
    logging.info(assmsg)


if __name__ == "__main__":
    main()
