cimport cython
@cython.boundscheck(False)

# def hamming_distance(seq1, seq2, limit=0):
#     i = 0
#     sum = 0
#     for i in range(len(seq1)):
#         if seq1[i] != seq2[i]:
#             sum += 1
#         if limit > 0 and sum > limit:
#             return limit + 1
#     return sum

cpdef int hamming_distance_cpdef(str seq1, str seq2, int limit=0):
    cdef int i, sum_
    i = 0
    sum_ = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            sum_ += 1
        if limit > 0 and sum_ > limit:
            return limit + 1
    return sum_

# cpdef int hamming_distance(a, b):
#     cdef char * aa = a
#     cdef char * bb = b
#     cdef int k, l, c
#     c = 0
#     l = len(a)
#     for k from 0 <= k < l:
#         if aa[k] != bb[k]:
#             c += 1
#     return c

# def get_kmers(sequence,k):
#     kmers = []
#     prev = ''
#     kmer = 0
#     offset = 0
#     width = len(sequence)
        
#     for i in range(width):
#         x = sequence[i]

#         if i==0 or x==prev: #if starting or staying on same kmer
#             if kmer < k:
#                 kmer+=1
#         else: #if ending a kmer
#             if kmer > 1:
#                 key = (kmer, prev, offset)
#                 if prev !='N': kmers.append(key)
#                 offset = i
#             kmer=1 
#         if i == width-1: # if end of the string
#             if kmer > 1:
#                 key = (kmer, x, offset)
#                 if x !='N': kmers.append(key)
#                 offset = i
#         prev=x
#     return kmers

# kmers_cython.pyx
# cpdef list get_kmers_cpdef(str sequence, int k):
#     cdef list kmers = []
#     cdef int prev = ''
#     cdef int kmer = 0
#     cdef int offset = 0
#     cdef int width = len(sequence)
#     cdef int x

#     cdef int i = 0
    
#     for i in range(width):
#         print(i)
#         x = sequence[i]

#         if i == 0 or x == prev:  # if starting or staying on the same kmer
#             if kmer < k:
#                 kmer += 1
#         else:  # if ending a kmer
#             if kmer > 1:
#                 key = (kmer, prev, offset)
#                 if prev != 'N':
#                     kmers.append(key)
#                 offset = i
#             kmer = 1
#         if i == width - 1:  # if end of the string
#             if kmer > 1:
#                 key = (kmer, x, offset)
#                 if x != 'N':
#                     kmers.append(key)
#                 offset = i
#         prev = x

#     return kmers

# kmers_cython.pyx
cpdef list get_kmers_cpdef(str sequence, int k):
    cdef list kmers = []
    cdef str prev = ""
    cdef int kmer = 0
    cdef int offset = 0
    cdef int width = len(sequence)

    cdef int i = 0
    
    cdef str x

    for i in range(width):
        x = sequence[i]
        # print(x)
        # print(prev)

        
        if i == 0 or x == prev:  # if starting or staying on the same kmer
            if kmer < k:
                kmer += 1
        else:  # if ending a kmer
            if kmer > 1:
                key = (kmer, prev, offset)
                if prev != "N":
                    kmers.append(key)
                offset = i
            kmer = 1
        if i == width - 1:  # if end of the string
            if kmer > 1:
                key = (kmer, x, offset)
                if x != "N":
                    kmers.append(key)
                offset = i
        prev = x

    return kmers


# def get_candidates(sequence, kmers_dict):
#     k=2
#     sequence_kmers = get_kmers(sequence,k)
#     pool = {}
#     for key in sequence_kmers:
#         kmer_len = key[0]
#         nt = key[1]
#         pos = key[2]
        
#         if key in kmers_dict:
#             pool.update({key: kmers_dict[key]})
        
#         up_pos = pos
#         down_pos = pos
        
        
#         while up_pos < k + pos:
#             new_key = (kmer_len,nt,up_pos)
#             if new_key in kmers_dict and new_key not in pool:
#                 pool.update({new_key: kmers_dict[new_key]})
#             new_key = (kmer_len+1,nt,up_pos)
#             if new_key in kmers_dict and new_key not in pool:
#                 pool.update({new_key: kmers_dict[new_key]})
#             up_pos+=1

#         thresh = max(0, down_pos-k)
#         while down_pos >= thresh:
#             new_key = (kmer_len,nt,down_pos)
#             if new_key in kmers_dict and new_key not in pool:
#                 pool.update({new_key: kmers_dict[new_key]})
#             new_key = (kmer_len+1,nt,down_pos)
#             if new_key in kmers_dict and new_key not in pool:
#                 pool.update({new_key: kmers_dict[new_key]})
#             down_pos-=1
        
#         unique_barcode_set = set()
#         for lst in pool.values():
#             unique_barcode_set.update(lst)
#         # end = time.time()
#         # print(end - start)
#         #print(len(unique_barcode_set))
#     return list(unique_barcode_set)

# candidates_cython.pyx
from libc.stdlib cimport malloc

cpdef list get_candidates_cpdef(str sequence, dict kmers_dict):
    cdef int k = 2
    cdef list sequence_kmers = get_kmers_cpdef(sequence, k)
    cdef dict pool = {}

    cdef int kmer_len
    cdef str nt
    cdef int pos

    cdef int up_pos 
    cdef int down_pos

    cdef tuple new_key

    cdef int thresh

    for key in sequence_kmers:
        kmer_len = key[0]
        nt = key[1]
        pos = key[2]

        if key in kmers_dict:
            pool[key] = kmers_dict[key]

        up_pos = pos
        down_pos = pos

        while up_pos < k + pos:
            new_key = (kmer_len, nt, up_pos)
            if new_key in kmers_dict and new_key not in pool:
                pool[new_key] = kmers_dict[new_key]
            new_key = (kmer_len + 1, nt, up_pos)
            if new_key in kmers_dict and new_key not in pool:
                pool[new_key] = kmers_dict[new_key]
            up_pos += 1

        thresh = max(0, down_pos - k)
        while down_pos >= thresh:
            new_key = (kmer_len, nt, down_pos)
            if new_key in kmers_dict and new_key not in pool:
                pool[new_key] = kmers_dict[new_key]
            new_key = (kmer_len + 1, nt, down_pos)
            if new_key in kmers_dict and new_key not in pool:
                pool[new_key] = kmers_dict[new_key]
            down_pos -= 1

    cdef set unique_barcode_set = set()
    for lst in pool.values():
        unique_barcode_set.update(lst)

    # Create Python strings directly without manual memory management
    cdef list result_list = []

    cdef bytes barcode_bytes 
    cdef char* c_barcode 
    cdef str newcode

    print(unique_barcode_set)
    
    for barcode in unique_barcode_set:
        barcode_bytes = barcode.encode('utf-8')
        c_barcode = <char*>malloc(len(barcode_bytes) + 1)
        c_barcode[:len(barcode_bytes)] = memoryview(barcode_bytes)
        c_barcode[len(barcode_bytes)] = '\0'
        #newcode = str(PyUnicode_FromStringAndSize(c_barcode, len(barcode_bytes)))
        newcode = str(<bytes>c_barcode)
        #print(newcode)
        result_list.append(newcode)

    return result_list



# def score_barcode_for_dict(seq, barcodes, max_edit_distance):
#     """
# 	this function scores a given sequence against all the barcodes. It returns the winner with Ns included.
# 	"""

#     if seq in barcodes:  # no need to check all barcodes
#         winner = seq  # barcode WITH Ns included
#     # elif min_score == len(barcodes_no_N):  # i.e. no matches allowed, and seq not in barcodes
#     #     winner = "no_match"
#         # print(seq)
#         # print(winner)
#         # print("")
#     else:  # mismatches allowed so need to check
#         dists = {}

#         for this_bc in barcodes:
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
#             #print(dists)
#             winner = filtered[0]  # barcode WITH Ns included
#     #     print(min_dist)
#     # print(seq)
#     # print(winner)
#     # print("")
#     return winner