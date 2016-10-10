import numpy as np
import pickle

SU_KEYWORDS = [
	('tracl', np.int32),
	('tracr', np.int32),
	('fldr', np.int32),
	('tracf', np.int32),
	('ep', np.int32),
	('cdp', np.int32),
	('cdpt', np.int32),
	('trid', np.int16),
	('nvs', np.int16),
	('nhs', np.int16),
	('duse', np.int16),
	('offset', np.int32),
	('gelev', np.int32),
	('selev', np.int32),
	('sdepth', np.int32),
	('gdel', np.int32),
	('sdel', np.int32),
	('swdep', np.int32),
	('gwdep', np.int32),
	('scalel', np.int16),
	('scalco', np.int16),
	('sx', np.int32),
	('sy', np.int32),
	('gx', np.int32),
	('gy', np.int32),
	('counit', np.int16),
	('wevel', np.int16),
	('swevel', np.int16),
	('sut', np.int16),
	('gut', np.int16),
	('sstat', np.int16),
	('gstat', np.int16),
	('tstat', np.int16),
	('laga', np.int16),
	('lagb', np.int16),
	('delrt', np.int16),
	('muts', np.int16),
	('mute', np.int16),
	('ns', np.uint16),
	('dt', np.uint16),
	('gain', np.int16),
	('igc', np.int16),
	('igi', np.int16),
	('corr', np.int16),
	('sfs', np.int16),
	('sfe', np.int16),
	('slen', np.int16),
	('styp', np.int16),
	('stas', np.int16),
	('stae', np.int16),
	('tatyp', np.int16),
	('afilf', np.int16),
	('afils', np.int16),
	('nofilf', np.int16),
	('nofils', np.int16),
	('lcf', np.int16),
	('hcf', np.int16),
	('lcs', np.int16),
	('hcs', np.int16),
	('year', np.int16),
	('day', np.int16),
	('hour', np.int16),
	('minute', np.int16),
	('sec', np.int16),
	('timbas', np.int16),
	('trwf', np.int16),
	('grnors', np.int16),
	('grnofr', np.int16),
	('grnlof', np.int16),
	('gaps', np.int16),
	('otrav', np.int16),
	## su variation
	('d1',np.float32),
	('f1',np.float32),
	('d2',np.float32),
	('f2',np.float32),
	('ungpow',np.float32),
	('unscale',np.float32),
	('ntr',np.int32),
	('mark',np.int16),
	('shortpad',np.int16),
	('unass',(np.int16,14)),
	## SEGY rev1
	#('cdpx', np.int32),
	#('cdpy', np.int32),
	#('Inline3D', np.int32),
	#('Crossline3D', np.int32),
	#('ShotPoint', np.int32),
	#('ShotPointScalar', np.int16),
	#('TraceValueMeasurementUnit', np.int16),
	#('TransductionConstantMantissa', np.int32),
	#('TransductionConstantPower', np.int16),
	#('TransductionUnit', np.int16),
	#('TraceIdentifier', np.int16),
	#('ScalarTraceHeader', np.int16),
	#('SourceType', np.int16),
	#('SourceEnergyDirectionMantissa', np.int32),
	#('SourceEnergyDirectionExponent', np.int16),
	#('SourceMeasurementMantissa', np.int32),
	#('SourceMeasurementExponent', np.int16),
	#('SourceMeasurementUnit', np.int16),
	#('UnassignedInt1', np.int32),
	#('UnassignedInt2', np.int32),
	]
SU_HEADER_DTYPE = np.dtype(SU_KEYWORDS)
SU_KEY_LIST = [keytype[0] for keytype in SU_KEYWORDS]

class SeismicTrace:
	def __init__(self,header,data,logs=[],nmo_picks={}):
		ntr,ns=data.shape
		# do not use data with only 1 trace
		hntr=header.shape[0]
		try:
			hns =header[0]['ns']
		except:
			hns =header['ns']
		if hntr != ntr:
			print("ntr of header and data do not match")
			return
		if hns != ns:
			print("ns of header and data do not match")
			return
		self.header = header.copy()
		self.data   = data.copy()
		self.log    = list(logs)
		self.nmo_picks  = nmo_picks

	def add_log(self,msg):
		self.log.append(msg)
	
	def print_log(self,nmo=False):
		for log in self.log:
			print(log)
		if nmo:
			print("nmo picks")
			print(self.nmo_picks)
	
	def logs(self):
		return list(self.log)

	def __sub__(self,trc):
		# substract trace (same size: ntr,ns)
		# output = current - trc
		out=SeismicTrace(self.header,self.data-trc.data,self.logs())
		out.add_log('sub')
		return out

	def __add__(self,trc):
		# add trace (same size: ntr,ns)
		# output = current + trc
		out=SeismicTrace(self.header,self.data+trc.data,self.logs())
		out.add_log('add')
		return out

	def write(self,filename):
		self.add_log("write: "+filename)
		with open(filename,'wb') as f:
			pickle.dump(self,f,pickle.HIGHEST_PROTOCOL)

	def write_su(self,filename):
		ntr,ns=self.data.shape
		su_file_dtype = np.dtype(SU_HEADER_DTYPE.descr + [('trace',('f4',ns))])
		output=np.empty(ntr,dtype=su_file_dtype)
		for key in SU_KEY_LIST:
			output[key]=self.header[key]
		output['trace']=self.data.astype(np.float32)
		output.tofile(filename)

def read(filename):
	with open(filename,'rb') as f:
		trace=pickle.load(f)
	trace.add_log("read: "+filename)
	return trace

def read_su(sufile):
	raw = open(sufile,'rb').read()
	# get ns
	su_header = np.fromstring(raw, dtype=SU_HEADER_DTYPE, count=1)
	ns = su_header['ns'][0]
	# define structure of a structured array
	su_file_dtype = np.dtype(SU_HEADER_DTYPE.descr + [('trace',('f4',ns))])
	# read traces
	su_traces = np.fromstring(raw, dtype=su_file_dtype)
	# return SeismicTrace
	header = su_traces[SU_KEY_LIST]
	data   = su_traces['trace'].astype(np.float64)
	out=SeismicTrace(header,data,['read_su: '+sufile])
	return out

def get_key(self,keyword):
	# return keyword values (with duplicated values)
	return self.header[keyword].copy()

def get_keys(self,keywords):
	# get header keywords
	# keywords: list/tuple of keywords to extract
	# output: value array
	out=[]
	for key in flatlist(keywords):
		out.append(get_key(self,key))
	return out

def get_key_unique(self,keyword):
	# return keyword values (without duplicated values)
	tmp=get_key(self,keyword)
	_, idx = np.unique(tmp, return_index=True)
	return tmp[np.sort(idx)]

def get_key_count(self,keyword):
	# return number of values for a keyword in SeismicTrace
	return len(get_key_unique(self,keyword))

def get_nshot(self,keyword="fldr"):
	# return number of shots in SeismicTrace
	# keyword="fldr": keyword to seperate shots
	return get_key_count(self,keyword)

def get_ns(self):
	# return ns (number of samples) of SeismicTrace
	try:
		ns=self.header[0]["ns"]
	except: # one trace data
		ns=self.header["ns"]
	return ns

def get_dt(self):
	# return dt of SeismicTrace [seconds]
	try:
		dt=self.header[0]["dt"]
	except: # one trace data
		dt=self.header["dt"]
	return dt*1.e-6

def get_ntr(self):
	# return number of traces in SeismicTrace
	if self.header.shape == ():
		return 0
	else:
		return self.header.shape[0]

def ntr_per_shot(self,keyword="fldr"):
	# return number of traces in each shot gather
	# keyword="fldr": keyword used to seperate each shot
	# output: array containing number of traces, size=nshot
	nshot=get_nshot(self,keyword)
	ntr=get_ntr(self)
	ntr_per_csg=np.zeros(nshot,dtype=np.int32)
	fldrs=get_key(self,keyword)
	table={}
	for ishot,fldr in zip(range(nshot),get_key_unique(self,keyword)):
		table[fldr]=ishot
	for itr in range(ntr):
		ishot = table[fldrs[itr]]
		ntr_per_csg[ishot] += 1
	return ntr_per_csg

def flatlist(val,dtype=False):
	# generate iterable from val (scalar or list)
	if dtype:
		return np.array([val],dtype=dtype).flatten()
	else:
		return np.array([val]).flatten()
	
def window(self,keyword,lst):
	# window SeismicTrace by keyword
	# keyword: keyword to window SeismicTrace
	# lst: list of keyword values to extract
	window = np.zeros(get_ntr(self),dtype=np.bool)
	key = get_key(self,keyword)
	for i in flatlist(lst,np.int32):
		window += key == i
	out=SeismicTrace(self.header[window], self.data[window],self.logs(),self.nmo_picks)
	out.add_log('window: key=%s range=%s'%(keyword,lst))
	return out

def print_range(self):
	# Print ranges of SeismicTrace keywords (non-zero values only)
	print("keyword:: min(loc) ~ max(loc): [first - last]// # of unique keys\n")
	for key in SU_KEY_LIST:
		keys=get_key(self,key)
		mn=keys.min()
		mx=keys.max()
		if mn==0 and mx==0:
			continue
		imn=keys.argmin()
		imx=keys.argmax()
		count=get_key_count(self,key)
		print("%s:: %d(%d) ~ %d(%d): [%d - %d]// %d"%(key,mn,imn,mx,imx,keys[0],keys[-1],count))

def trace_split(self,keyword):
	# split SeismicTrace by keyword (SeismicTrace doesn't need to be sorted)
	# keyword: keyword to split SeismicTrace
	out=[]
	vals=get_key_unique(self,keyword)
	for val in vals:
		out.append(window(self,keyword,val))
	return out

def trace_sort(self,keys):
	# Sort SeismicTrace by keywords
	# keys: sort keys (ex. ['+fldr','-offset'])
	#      + or no sign: ascending order
	#      - sign: descending order
	# output: sorted SeismicTrace
	sortkeys=list(keys)
	for i,key in enumerate(keys):
		order=1
		if key[0] in ['+','-']:
			order = int(key[0]+'1')
			k=key[1:]
		else:
			k=key
		sortkeys[i]=get_key(self,k)*order
	sortkeys.reverse()
	ind = np.lexsort(sortkeys)
	out=SeismicTrace(self.header[ind],self.data[ind],self.logs(),self.nmo_picks)
	out.add_log("sort: keys=%s"%(','.join(list(keys))))
	return out

