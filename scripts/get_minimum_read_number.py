if __name__ == '__main__':

    import argparse, os
    parser = argparse.ArgumentParser(description='Gets minimum number from the last column from input txt and returns it')
    parser.add_argument('--input', dest='input_file', type=str, help='input_file')
    args = parser.parse_args()
    f = open(args.input_file,'r')
    data = []
    for l in f.readlines():
        data.append(int(l.split()[-1]))
    minimum = min(data)
    print(minimum)	 
