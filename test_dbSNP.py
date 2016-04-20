#!/usr/bin/python
import argparse
import os
import re
import pdb
#from guppy import hpy
#import resource
#from profilestats import profile
#import cProfile
#import line_profiler
#import time


def vcf2json_alt(v, conv,vcf_content ):
	vf = open(v, 'r')
	status_bar = ''
	data_pool = 10000
	
	for aline in vf:
		if(not re.match('^#', aline)):
			aline = re.sub('\s+$', '', aline)
			content = aline.split('\t')
			aVar = {}
			full_loc = conv[content[0]]
			alt_loc = conv[content[0]].split('_')
			var_chr = alt_loc[0]
			
			var_ncbi = ''
			if len(alt_loc) >1:
				var_ncbi = alt_loc[1]
			var_ct = ''
			if len(alt_loc) >2:
				var_ct = alt_loc[2]
			var_p   = content[1]
			#var_rsid  = content[2]
			var_ref = content[3]
			var_alt = content[4]
			if( var_chr == 'chrX' or var_chr == 'X'):
				var_chr = 'chr23'
			elif (var_chr =='chrY' or var_chr == 'Y'):
				var_chr = 'chr24'
			elif(var_chr =='MT'):
				var_chr = 'chr25'
			if status_bar != full_loc:
				status_bar = full_loc
				print(status_bar)
				#print(status_bar, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, ' MB')
			if var_ct == '':
				var_id = var_chr+'_'+var_ncbi+':'+var_p
			else:
				var_id = var_chr+'_'+var_ncbi+'_'+var_ct+':'+var_p
			alt_allele = re.split(',', var_alt)
			num_allele = len(re.split(',', var_alt))
			grp_size = []
			total_size = 0
			cln_allele = {}
			info = content[-1].split(';')
			
			#vcf_content = {}
			tmp = {}
			for ele in info:  ## parse the info column of VCF file
				pair = ele.split('=')
				
				if(len(pair) ==2):
					tmp[pair[0]] = pair[1]
			tmp['o'] = var_alt
			tmp['r'] = var_ref
			#if 'rs'+tmp['RS'] in ucsc.keys():
			#	tmp['ucsc'] = 1
			#else:
			#	tmp['ucsc'] = 0
			
			if var_id not in vcf_content.keys():
				vcf_content[var_id] = [tmp]
			else:
				vcf_content[var_id].append(tmp)
	vf.close()
	return(status_bar)

def inUCSC(fl, ucsc):
	uf = open(fl, 'r')
	header = uf.readline()
	for u in uf:
		aline = u.strip().split('\t')
		ucsc[aline[3]] = 1

#@profile
#@profile
def random_sort2(n):
	l = [random.random() for i in range(n)]
	l.sort()
	return l



ucsc = {}
inUCSC('GRCh38_ucsc_commonSNPs144.bed', ucsc)


def testdbSNP(infle, ucsc):
	conv ={}
	convt = open('GCF_000001405.28.assembly.txt', 'r')
	for cl in convt:
		if not cl.startswith('#'):
			aline = cl.strip().split('\t')
			conv[aline[6]] = aline[-1]
	convt.close()

	#infle = args[1]

	oridata = {}
	vcf_small_file = vcf2json_alt(infle, conv, oridata)
	#vcf2json_alt('testvcf', conv, oridata)
	#out_f = open(infle+'_test.json', 'w')
	out_f = open(vcf_small_file.strip()+'_test.json', 'w')

	#meta = []
	for var,all_info in oridata.items():
		id_inf = ''
		id_json = ''
		#var_pos = ''
		delim = re.compile('[:_]')
		id_e = delim.split(var)
		id_e[0] = re.sub('chr', '',id_e[0])
		var_pos = id_e[-1]
		if id_e[0] == 'Un':
			id_e[0] = '27'
		if len(id_e) == 4:
			#id_e = id_json.split('_')
			#id_e = re.sub('chr', '',id_e)
			#var_pos = id_e[3]
			id_json = '{\"_id\":{\"c\":'+id_e[0]+',\"ct\":\"'+id_e[2]+'\",\"ncbi\":\"'+ id_e[1] + '\",\"p\":'+id_e[3]+'},'
		elif len(id_e) == 3:
			#var_pos = id_e[2]
			id_json = '{\"_id\":{\"c\":'+id_e[0]+',\"ncbi\":\"'+ id_e[1] + '\",\"p\":'+id_e[2]+'},'
		elif len(id_e) ==2:
			#var_pos = id_e[1]
			id_json = '{\"_id\":{\"c\":'+id_e[0]+',\"p\":'+id_e[1]+'},'

		#json_line = id_json +'\"f\":[{'
		json_line = []
		for info in all_info:
			meta = []
			meta.append('\"i\":'+'\"rs'+info['RS'] +'\"')
			alt_allele = info['o'].split(',')
			f_o = []
			for a in alt_allele:
				f_o.append('\"'+a+'\"')
			meta.append('\"o\":[' + ','.join(f_o) +']')
			if (info['RSPOS'] != var_pos and info['VC'] != 'DIV') or (info['VC'] == 'DIV' and int(info['RSPOS']) != (int(var_pos)+1)):
				meta.append('\"opos\":'+ info['RSPOS'])
			meta.append('\"r\":'+'\"'+info['r'] +'\"')
			meta.append('\"s\":'+info['SAO'])
			info['VC'] = re.sub('DIV', 'INDEL', info['VC'])
			meta.append('\"t\":['+'\"'+info['VC'] +'\"]')
			if 'rs'+info['RS'] in ucsc.keys():
			#if info['ucsc'] ==1:
				meta.append('\"ucsc\":1')

			#json_line = json_line + ',{' + ','.join(meta)+ '}'
			json_line.append('{' + ','.join(meta) + '}')
		out_line = id_json +'\"f\":[' + ','.join(json_line) + ']}'
		out_f.write(out_line+'\n')

	out_f.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input')
	args = parser.parse_args()
	infls = open(args.input, 'r')
	
	ucsc = {}
	inUCSC('GRCh38_ucsc_commonSNPs144.bed', ucsc)

	for f in infls:
		testdbSNP(f.strip(), ucsc)
	infls.close()
	ucsc = {}

