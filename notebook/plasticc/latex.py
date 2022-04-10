import glob


ll = open('tmp.txt', 'w')
for f in glob.glob('*png'):
    objid = f.replace('.png','')
    ll.write (
        '\\begin{figure}[H]\n'+\
        '\centering\n'+\
        '\includegraphics[width=\\textwidth]{fits/%s.png}\n'%objid+\
        '\caption{objid=%s}\n'%objid+\
        '\end{figure}\n\n'
    )
ll.close()
