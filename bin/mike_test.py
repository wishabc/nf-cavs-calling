

def find_min_and_max(blocks):
    min_block_ind, min_block_value = -1, 10**10
    max_block_ind, max_block_value = -1, -1
    for i, x in enumerate(blocks):
        if x < min_block_value:
            min_block_value = x
            min_block_ind = i
        if x > max_block_value:
            max_block_value = x
            max_block_ind = i
        
    return min_block_ind, max_block_ind


def get_minimum_num_moves(blocks):
    min_index, max_index = find_min_and_max(blocks)
    if min_index > max_index:
        # we are going to change the order of min and max with one swap,
        # -1 to not count twice
        n = min_index + (len(blocks) - 1 - max_index) - 1
    else:
        n = min_index + (len(blocks) - 1 - max_index)
    return n

if __name__ == '__main__': 
    # blocks = blocks = [3,2,1]
    # blocks = [4,11,9,10,12]
    # blocks = [2,4,3,1,6]
    # blocks = [2,4,3,6,1]
    # blocks = [1,4,3,6,2]
    # blocks = [10,4,30,6,2]
    blocks = [100,4,30,6,2]
    #blocks = [100, 4, 30, 6, 2, 1, 2, 10, 3]
    # 1) [2 3 1 ]
    # 2) [2,1,3]
    # 3) [1,2,3] return 3
    print(get_minimum_num_moves(blocks))
