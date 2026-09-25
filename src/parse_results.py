# Produce data for 2 parameter tikz plot 
import re, glob
import numpy as np


def is_valid_template(path):
    '''* for all'''
    subs = path.split('*')
    return all(len(sub) > 0 for sub in subs)


def get_param_value(path, name, pattern=None, which_num=0):
    '''name_L{*}'''
    numbers = re.compile(r'[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?')
    if pattern is None:
        pattern = re.compile(rf'\w+{name}[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?')
    else:
        pattern = re.compile(pattern)
    try: 
        matched, = pattern.findall(path)
    except ValueError:
        return None

    num = matched[matched.index(name)+len(name):]
    num = numbers.findall(num)[which_num]
    
    return float(num)


def get_data(path, col_name, normalize=False, skip=0, preprocess=None):
    '''Extract column'''
    with open(path) as lines:
        for k in range(skip):
            next(lines)
        header = next(lines).strip().split()
        try:
            col = header.index(col_name)-1
            data = [float(line.strip().split()[col]) for line in lines]
            if preprocess is not None:
                data = [preprocess(val) for val in data]
        except ValueError:
            if col_name == 'level':
                data = [l for (l, line) in enumerate(lines)]
            else:
                print(col_name, header, '<<<<')
                exit()
        if normalize:
            data = np.array(data)
            data /= data[0]
        
        return data

    
def contains(iterable, num, tol=1E-14):
    '''Is num in iterable?'''
    for i, item in enumerate(iterable):
        if abs(item-num) < tol:
            return i
    return None


def extract_data(template, variable_ranges, xcol, ycol, normalize, reverse_sort_x=False, skip=0, pattern=None, which_num=0, preprocess=None):
    '''Look for data in tables and align'''
    assert len(variable_ranges) == template.count('*')

    paths = dict()
    variables = tuple(variable_ranges.keys())
    print(template)
    print(glob.glob(template))
    for path in glob.glob(template):
        
        if 'tikz' in path: continue

        print('\t', path)
        key = ()
        is_valid = True
        for variable in variables:
            value = get_param_value(path, variable, pattern=pattern, which_num=which_num)
            if value is None:
                break
            variable_range = variable_ranges[variable]
            # print('->', variable, value, variable_range)
            is_valid = not variable_range or contains(variable_range, value) is not None
            key = key + (value, )
            if not is_valid:
                break
        # Only valid can asssign
        if is_valid:
            paths[key] = path

        print(path, key, variables, is_valid)
            
    # Want to align data for union of all indep
    X = set()
    for key in paths:
        path = paths[key]
        X.update(get_data(path, xcol, skip=skip, preprocess=preprocess))
    X = np.array(sorted(X, reverse=reverse_sort_x))

    aligned_y = dict()
    # Now we insert NaNs as indicator of missing data
    for key in paths:
        path = paths[key]
        x = get_data(path, xcol, skip=skip, preprocess=preprocess)
        # Make room for larger
        Y = np.nan*np.ones_like(X)
        y = get_data(path, ycol, normalize, skip=skip, preprocess=preprocess)
        for xi, yi in zip(x, y):
            # Look for where to insert
            idx = contains(X, xi)
            if idx is not None:
                Y[idx] = yi
        # Final data
        aligned_y[key] = Y

    return X, aligned_y, variables


def compute(directory, template, ranges, xcol, ycol, normalize, sort_x, skip, pattern, which_num):
    '''Do it'''
    template = os.path.join(directory, template)

    if 'cond' in ycol:
        preprocess = lambda v: round(v, 3) if v < 100 else float(int(v))
    else:
        preprocess = None
    
    X, aligned_y, variables = extract_data(template=template,
                                           variable_ranges=ranges,
                                           xcol=xcol,
                                           ycol=ycol,
                                           normalize=normalize,
                                           reverse_sort_x=sort_x,
                                           skip=skip,
                                           preprocess=preprocess,
                                           pattern=pattern,
                                           which_num=which_num)
    # assert aligned_y
    print('X', X)
    print(aligned_y, variables)
    # Now for printing in tikz
    data = np.column_stack([X] + [aligned_y[key] for key in aligned_y])

    def tikz_column_name(key, variables=variables):
        return ''.join([f'{var}{keyi:.4E}' for (var, keyi) in zip(variables, key)])
    
    header = ' '.join(['x'] + [tikz_column_name(key) for key in aligned_y])

    tikz_path, ext = os.path.splitext(template.replace('*', 'varied'))
    tikz_path = '_'.join([tikz_path, 'paramRobust',
                          f'x{args.x_column.upper()}', f'y{y_column.upper()}', f'NORMALIZED{args.normalize}'
                          'tikz'])
    tikz_path = f'{tikz_path}.txt'

    with open(tikz_path, 'w') as out:
        out.write('%s\n' % header)
        np.savetxt(out, data)

    # Self inspection
    prev = None
    for key in sorted(aligned_y):
        if key:
            print(key, '->', np.round((aligned_y[key])[np.isfinite(aligned_y[key])], 2))

    import tabulate
    headers = sorted(aligned_y)

    keys = sorted(tuple(aligned_y.keys()), reverse=False)
    table = np.vstack([key for key in keys])

    ultima = np.array([aligned_y[key][np.where(~np.isnan(aligned_y[key]))[0][-1]] for key in keys])
    disp_ultima = np.where(np.abs(ultima) > 50, np.round(ultima,0), np.round(ultima, 2))

    if all(len(aligned_y[key]) > 3 for key in keys):
        penultima = np.array([aligned_y[key][np.where(~np.isnan(aligned_y[key]))[0][-2]] for key in keys])
        rel = np.abs(ultima-penultima)/ultima

        table = np.c_[table, ultima, np.abs(ultima-penultima), rel]

        headers = variables + ('data', 'diff', 'rel')
        print()
        print(tabulate.tabulate(table, headers=headers))
        print()
    else:
        table = np.c_[table, disp_ultima]

        headers = variables + ('data', )
        print('xxxxxxxxxx')
        print(tabulate.tabulate(table, headers=headers))
        print()            

    table = np.vstack([key for key in keys])
    ultima = np.array([aligned_y[key][np.where(np.logical_or(~np.isnan(aligned_y[key]),
                                                             np.isnan(aligned_y[key])))[0][:]] for key in keys])
    
    table = np.c_[table, ultima]
    nrows, ncols = table.shape
    for row in table:
        try:
            fmt = ' & '.join(['%g'] + ['%2.2f']*(ncols-1)) + r'\\'
            print(fmt%tuple(row))
        except:
            row = [row[0]] + row[1].tolist()
            fmt = ' & '.join(['%g'] + ['%2.2f']*(len(row)-1)) + r'\\'
            print(fmt%tuple(row))
    print
    
    print(os.path.abspath(tikz_path))

    return table


def value_fmt(y):
    return {'niters' : '%d',
            'level' : '%d',
            'condKSP': '%2.1f'}.get(y, '%g')

# --------------------------------------------------------------------

if __name__ == '__main__':
    import argparse, os
    import numpy as np
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Where to get data
    parser.add_argument('directory', type=str, help='Directory with files')

    parser.add_argument('-x_column', type=str, default='ndofs', help='Name of independent variable')
    parser.add_argument('-y_column', type=str, nargs='+', help='Name of dependent variable')        

    parser.add_argument('-normalize', type=int, default=0, choices=(0, 1))


    parser.add_argument('-precond', type=str, default='standard')
    parser.add_argument('-inverseQ', type=str, default='lu')
    
    # bcs
    parser.add_argument('-ubcs', type=str, help='Spec of bcs for momentum: D or T', default='TTDD')
    parser.add_argument('-pbcs', type=str, help='Spec of bcs for mass: P or F', default='FFFF')    
    # Material
    parser.add_argument('-alpha', type=float, default=1E0)
    parser.add_argument('-K', type=float, default=1)
    parser.add_argument('-mu', type=float, default=0.5)
    parser.add_argument('-lmbda', type=float, default=1)
    parser.add_argument('-c', type=float, default=0.0)    
    
    args, _ = parser.parse_known_args()

    template = f'precond{args.precond}/cvrg_L*_PRECOND{args.precond}_INVERSEQ{args.inverseQ}_UBCS{args.ubcs}_PBCS{args.pbcs}_ALPHA{args.alpha}_K{args.K}_MU{args.mu}_LMBDA{args.lmbda}_C{args.c}.txt'

    pattern = None
    which_num = 0

    variable_ranges = {#'LENGTH': (1, 2, 5, 10, 20, 50, 100)}
        'L': (1, 5, 10, 20, 40, 80)
    }    
    
    reverse_sort_x = {'hmin': True, 'h': True, '|ep|_0': True,
                      'ndofs': False,
                      'level': False}[args.x_column]

    tables = []
    for y_column in args.y_column:
        table = compute(directory=args.directory,
                        template=template,
                        ranges=variable_ranges,
                        xcol=args.x_column,
                        ycol=y_column,
                        normalize=args.normalize,
                        sort_x=reverse_sort_x,
                        skip=0,
                        pattern=pattern,
                        which_num=which_num)
        tables.append(table)
    (nrows, ncols), = set(table.shape for table in tables)
    ntables = len(tables)

    for r in range(nrows):
        value = value_fmt(args.y_column[0])
        if ntables > 1:
            value = value + '(' + ', '.join(map(value_fmt, args.y_column[1:])) + ')'
        fmt = ' & '.join(['%g'] + [value]*(ncols-1)) + r'\\'
        row_ = (tables[0][r, 0], ) + sum(zip(*(t[r, 1:] for t in tables)), ())
        row_ = tuple(v if not np.isnan(v) else -1 for v in row_)
        print(fmt % tuple(row_))
    print
