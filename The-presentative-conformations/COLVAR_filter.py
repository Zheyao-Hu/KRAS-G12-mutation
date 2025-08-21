import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Filter PDB file based on user criteria.")
    parser.add_argument("-i", "--input", "--colvar", type=str, default='COLVAR', help="Input COLVAR file, default=COLVAR")
    parser.add_argument("-o", '--output', type=str, default="COLVER_filtered", help='Output file, default=COLVAR_filtered')
    parser.add_argument("--set", nargs="*", 
            help="Custom params, such as:\na=2,3 means 2<=a<3;\nb=oo,6 means b<=6;\nc=-2,oo means c>=-2"
    )
    return parser.parse_args()

def is_right(val):
    if type(val) == str and val.count(',')==1:
        return True
    return False

def parse_value(val):
    con1 = ',' in val
    con2 = 'oo,' in val
    con3 = ',oo' in val

    if con1:
        pa = val.split(',')
        if not con2 and not con3:
            return [float(pa[0]), float(pa[1])]
        if con2 and not con3:
            return [-float('inf'), float(pa[1])]
        if not con2 and con3:
            return [float(pa[0]), float('inf')]
    return [-float('inf'), float('inf')]

def iter_file(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            yield line.strip('\n')

def main():
    args = parse_args()
    user_vars = {}
    if args.set:
        for item in args.set:
            if '=' in item:
                k, v = item.split('=', 1)
                if is_right(v):
                    user_vars[k] = parse_value(v)

    print(f'- Custom paras:{user_vars}')

    def judge_value(title, line_float_list):
        for k,v in user_vars.items():
            if not v[0] <= line_float_list[title.index(k)] <= v[1]:
                return False
        return True
    
    title = []
    i, j = 0, 0
    with open(args.output, 'w', encoding='utf-8') as fw:
        for line in iter_file(args.input):

            if line.startswith('#!'):
                if 'FIELDS' in line:
                    linel = line.split()
                    title.extend(linel[2:])
                    print(f"- title: {title}")
                print(line, file=fw)
                continue
            
            i += 1
            linell = list(map(float, line.split()))
            if judge_value(title, linell):
                print(line, file=fw)
                j += 1

    print(f"- Title index {[title.index(k) for k in user_vars.keys()]} was/were selected")
    print(f'- There are {i} data lines in total, of which {j} are retained.')

if __name__ == "__main__":
    main()
    print('Done')

