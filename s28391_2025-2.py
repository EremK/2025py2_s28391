from Bio import Entrez as E,SeqIO as S
import matplotlib.pyplot as p,sys,time
q=input
E.email=q('mail: ');k=q('key: ');E.api_key=k;E.tool='x'
tx=q('taxid: ').strip()or sys.exit()
mn=int(q('min: ')or 0);mx=int(q('max: ')or 1e9)

while 1:
  try:
    j=E.read(E.esearch(db='nucleotide',
        term=f'txid{tx}[Organism]',usehistory='y'));break
  except:print('.',end='',flush=True);time.sleep(1)

n=int(j['Count'])or sys.exit('no rec')
w,j=j['WebEnv'],j['QueryKey']
bag=[];t=time.time()
for s in range(0,n,500):
  with E.efetch(db='nucleotide',rettype='gb',retmode='text',
      retstart=s,retmax=500,webenv=w,query_key=j)as f:
    for r in S.parse(f,'genbank'):
      L=len(r.seq)
      if mn<=L<=mx:
        bag.append((r.id.split('.')[0],L,r.description.replace(',',' ')[:80]))
  if not k and time.time()-t<.34:time.sleep(.34)
bag or sys.exit('-')
bag.sort(key=lambda x:-x[1])
csv=f'r{tx}.csv'
with open(csv,'w')as o:
  o.write('acc,length,desc\n')
  o.writelines(f'{a},{l},{d}\n'for a,l,d in bag)

p.plot([a for a,_,_ in bag],[l for _,l,_ in bag],'o')
p.xticks(rotation=90,fontsize=6);p.tight_layout()
png=f'g{tx}.png';p.savefig(png,dpi=120)