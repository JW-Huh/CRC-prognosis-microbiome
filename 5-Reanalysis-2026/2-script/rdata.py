"""Minimal reader for R save()/save.image() RDX2/RDX3 XDR files -> dict of top-level objects (data.frames as pandas)."""
import gzip, struct, sys
import numpy as np, pandas as pd

NIL=254; REF=255; LISTSXP=2; SYMSXP=1; CHARSXP=9; STRSXP=16; INTSXP=13; REALSXP=14; LGLSXP=10; CPLXSXP=15; RAWSXP=24
VECSXP=19; EXPRSXP=20; CLOSXP=3; ENVSXP=4; PROMSXP=5; LANGSXP=6; SPECIALSXP=7; BUILTINSXP=8; DOTSXP=17; BCODESXP=21; EXTPTRSXP=22; WEAKREFSXP=23; S4SXP=25
EMPTYENV=242; BASEENV=241; GLOBALENV=253; UNBOUND=252; MISSINGARG=251; BASENAMESPACE=250; NAMESPACESXP=249; PACKAGESXP=248; PERSISTSXP=247
ALTREP=238; ATTRLISTSXP=239; ATTRLANGSXP=240; BCREPDEF=244; BCREPREF=243; CLASSREFSXP=246; GENERICREFSXP=245; EMPTYENV_=242

class RObj:
    __slots__=("type","value","attrs","tag")
    def __init__(s,t,v=None,a=None,tag=None): s.type=t; s.value=v; s.attrs=a; s.tag=tag
    def __repr__(s): return f"RObj({s.type})"

class Reader:
    def __init__(self, data, enc="cp949"):
        self.d=data; self.p=0; self.refs=[]; self.enc=enc; self.bcreps=None
    def int(self):
        v=struct.unpack(">i", self.d[self.p:self.p+4])[0]; self.p+=4; return v
    def ints(self,n):
        a=np.frombuffer(self.d, dtype=">i4", count=n, offset=self.p); self.p+=4*n; return a.astype(np.int64)
    def dbls(self,n):
        a=np.frombuffer(self.d, dtype=">f8", count=n, offset=self.p); self.p+=8*n; return a.astype(float)
    def charsxp(self, flags):
        n=self.int()
        if n==-1: return None
        b=self.d[self.p:self.p+n]; self.p+=n
        lv=flags>>12
        enc = "utf-8" if lv & 8 else ("latin-1" if lv & 4 else self.enc)
        try: return b.decode(enc)
        except: return b.decode("utf-8","replace")
    def item(self):
        flags=self.int(); t=flags&0xFF; has_attr=bool(flags&(1<<9)); has_tag=bool(flags&(1<<10)); lv=flags>>12
        if t==NIL: return None
        if t in (EMPTYENV,BASEENV,GLOBALENV,UNBOUND,MISSINGARG,BASENAMESPACE): return RObj(t)
        if t==REF:
            idx=flags>>8
            if idx==0: idx=self.int()
            return self.refs[idx-1]
        if t in (PERSISTSXP, PACKAGESXP, NAMESPACESXP):
            assert self.int()==0; n=self.int(); s=[self.item() for _ in range(n)]
            o=RObj(t,s); self.refs.append(o); return o
        if t==SYMSXP:
            s=self.item(); o=RObj(SYMSXP, s); self.refs.append(o); return o
        if t==ENVSXP:
            locked=self.int(); o=RObj(ENVSXP); self.refs.append(o)
            o.value=dict(enclos=self.item(), frame=self.item(), hashtab=self.item(), attrs=self.item()); return o
        if t in (LISTSXP,LANGSXP,CLOSXP,PROMSXP,DOTSXP,ATTRLISTSXP,ATTRLANGSXP):
            a=self.item() if has_attr else None
            tag=self.item() if has_tag else None
            car=self.item(); cdr=self.item()
            return RObj(t,(car,cdr),a,tag)
        if t==EXTPTRSXP:
            o=RObj(t); self.refs.append(o); o.value=(self.item(), self.item())
            if has_attr: o.attrs=self.item()
            return o
        if t==WEAKREFSXP:
            o=RObj(t); self.refs.append(o)
            if has_attr: o.attrs=self.item()
            return o
        if t in (SPECIALSXP,BUILTINSXP):
            n=self.int(); s=self.d[self.p:self.p+n].decode(); self.p+=n; return RObj(t,s)
        if t==CHARSXP: return self.charsxp(flags)
        if t in (LGLSXP,INTSXP):
            n=self.int(); v=self.ints(n)
            if t==LGLSXP: v=np.where(v==-2147483648, np.nan, v)
            o=RObj(t,v)
        elif t==REALSXP:
            n=self.int(); o=RObj(t,self.dbls(n))
        elif t==CPLXSXP:
            n=self.int(); v=self.dbls(2*n); o=RObj(t,v[0::2]+1j*v[1::2])
        elif t==STRSXP:
            n=self.int(); vals=[]
            for _ in range(n):
                f=self.int(); assert (f&0xFF)==CHARSXP, f
                vals.append(self.charsxp(f))
            o=RObj(t,vals)
        elif t in (VECSXP,EXPRSXP):
            n=self.int(); o=RObj(t,[self.item() for _ in range(n)])
        elif t==RAWSXP:
            n=self.int(); o=RObj(t,self.d[self.p:self.p+n]); self.p+=n
        elif t==S4SXP: o=RObj(t)
        elif t==BCODESXP:
            n=self.int(); saved=self.bcreps; self.bcreps=[None]*n; o=RObj(t,self.bc1()); self.bcreps=saved
        elif t==ALTREP:
            info=self.item(); state=self.item(); attr=self.item(); o=RObj(t,(info,state)); o.attrs=attr
            return self.expand_altrep(o)
        else:
            raise ValueError(f"unhandled type {t} at {self.p}")
        if has_attr: o.attrs=self.item()
        return o
    def bc1(self):
        code=self.item(); n=self.int(); consts=[]
        for _ in range(n):
            t=self.int()
            if t==BCODESXP: consts.append(self.bc1())
            elif t in (LANGSXP,LISTSXP,BCREPDEF,BCREPREF,ATTRLANGSXP,ATTRLISTSXP): consts.append(self.bclang(t))
            else: consts.append(self.item())
        return (code,consts)
    def bclang(self,t):
        if t==BCREPREF: return self.bcreps[self.int()]
        if t in (BCREPDEF,LANGSXP,LISTSXP,ATTRLANGSXP,ATTRLISTSXP):
            pos=-1
            if t==BCREPDEF: pos=self.int(); t=self.int()
            attr=self.item() if t in (ATTRLANGSXP,ATTRLISTSXP) else None
            tag=self.item(); car=self.bclang(self.int()); cdr=self.bclang(self.int())
            o=RObj(t,(car,cdr),attr,tag)
            if pos>=0: self.bcreps[pos]=o
            return o
        return self.item()
    def expand_altrep(self,o):
        info,state=o.value
        cls = info.value[0].value if info and info.type==LISTSXP else None
        name = cls.value if isinstance(cls,RObj) else cls
        if name=="compact_intseq":
            n,start,step=state.value; v=np.arange(start, start+step*n, step).astype(np.int64); r=RObj(INTSXP,v); r.attrs=o.attrs; return r
        if name=="compact_realseq":
            n,start,step=state.value; r=RObj(REALSXP,np.arange(start,start+step*n,step).astype(float)); r.attrs=o.attrs; return r
        if name in ("wrap_integer","wrap_real","wrap_string","wrap_logical"):
            r=state.value[0]; r.attrs=o.attrs if o.attrs is not None else r.attrs; return r
        if name=="deferred_string":
            r=state.value[0]; r.attrs=o.attrs; return r
        return o

def attrs_dict(o):
    d={}; a=o.attrs if isinstance(o,RObj) else None
    while isinstance(a,RObj) and a.type in (LISTSXP,ATTRLISTSXP):
        tag=a.tag.value if isinstance(a.tag,RObj) else a.tag
        d[tag]=a.value[0]; a=a.value[1]
    return d

def to_py(o):
    if o is None: return None
    if not isinstance(o,RObj): return o
    at=attrs_dict(o)
    if o.type in (INTSXP,REALSXP,LGLSXP,STRSXP,CPLXSXP):
        v=o.value
        if o.type==INTSXP: v=np.where(v==-2147483648, np.nan, v.astype(float)) if (v==-2147483648).any() else v
        if "levels" in at and o.type==INTSXP:   # factor
            lev=at["levels"].value; v=pd.Categorical.from_codes(np.nan_to_num(v,nan=0).astype(int)-1, lev) if not np.isnan(v).any() else pd.Series(v).map(lambda i: lev[int(i)-1] if not np.isnan(i) else None).values
        if "names" in at: return pd.Series(v, index=at["names"].value)
        return np.asarray(v) if not isinstance(v,(pd.Categorical,)) else v
    if o.type==VECSXP:
        cls=at.get("class"); cls=cls.value if isinstance(cls,RObj) else []
        if "data.frame" in cls:
            cols={n: to_py(c) for n,c in zip(at["names"].value, o.value)}
            df=pd.DataFrame({k:(v.values if isinstance(v,pd.Series) else v) for k,v in cols.items()})
            rn=at.get("row.names")
            if isinstance(rn,RObj) and rn.type==STRSXP: df.index=rn.value
            return df
        if "names" in at: return {n: to_py(c) for n,c in zip(at["names"].value, o.value)}
        return [to_py(c) for c in o.value]
    return o

def load(path):
    raw=open(path,"rb").read(); data=gzip.decompress(raw) if raw[:2]==b"\x1f\x8b" else raw
    assert data[:5] in (b"RDX2\n",b"RDX3\n"); r=Reader(data[5:])
    fmt=r.d[r.p:r.p+2]; r.p+=2
    ver=r.int(); wv=r.int(); mv=r.int()
    if ver==3:
        n=r.int(); r.enc=r.d[r.p:r.p+n].decode(); r.p+=n
    top=r.item(); out={}
    while isinstance(top,RObj) and top.type==LISTSXP:
        out[top.tag.value]=top.value[0]; top=top.value[1]
    return out

if __name__=="__main__":
    objs=load(sys.argv[1])
    for k,v in objs.items():
        print(k, getattr(v,"type",type(v)))
