import sys
import calculate_prevalences as cap
import os
import pandas as pd
import plotly.express as px
import warnings
warnings.filterwarnings('ignore')
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import yaml
import json
import plotly.graph_objs as go
import plotly.io as pio
import numpy as np
import PCA
import math
import subprocess

def get_mutation_counts(mutation_count_file):
	mutation_counts = pd.read_csv(mutation_count_file,
								  header=list(range(6)),
								  index_col=0)
	return mutation_counts

def get_mutation_coverage(mutation_coverage_file):
	mutation_coverage = pd.read_csv(mutation_coverage_file,
									index_col=0,
									header=list(range(6)))
	return mutation_coverage

def filter_mutations(mutation_count_file):
	mutation_counts = get_mutation_counts(mutation_count_file)
	all_columns=mutation_counts.keys()
	mutation_dict={column[2]:column_number for column_number, column in enumerate(all_columns) if column[3]=='missense_variant'}
	filtered_mutations=list(mutation_dict.keys())
	return filtered_mutations

def get_metadata_columns(metadata_table):
	metadata = pd.read_csv(metadata_table, sep='\t')
	columns = list(metadata.keys())
	return columns

def generate_sample_summary_dropdowns(metadata_table):
	metadata_columns = get_metadata_columns(metadata_table)
	sample_column = widgets.Dropdown(
		options=metadata_columns,
		# value='Sample_ID',
		description='Sample:',
		disabled=False,
	)
	summary_column = widgets.Dropdown(
		options=metadata_columns,
		description='Summary:',
		disabled=False,
	)
	return sample_column, summary_column

def generate_country_dropdown():
	country_shortcuts={
		'DRC': {
			"lat": -1.6815695315287824,
			"lon": 22.744896416745945,
			'zoom': 3.8
		},
		'Uganda': {
			"lat": 1.5,
			'lon': 32,
			'zoom': 5.6
		},
		'Tanzania': {
			"lat": -5.7,
			"lon": 35,
			'zoom': 5.2
		},
	}
	countries = list(country_shortcuts.keys())
	countries.append('custom')
	country = widgets.Dropdown(
		options=countries,
		description='Country:',
		disabled=False,
	)
	return country, country_shortcuts


def special_sort(mutations):
	'''
	Intersects the variants of interest with observed mutations from the AA tables
	default sort will sort mutations by gene and then reference amino acid name (e.g. alanine before lysine) - this sorter instead
	sorts by gene and then amino acid position
	'''
	sorting_list=[]
	for mutation in mutations:
		gene='-'.join(mutation.split('-')[:-1]).lower()
		pos=int(mutation.split('-')[-1][3:-3]) #strip out amino acid names, leaving only positions
		sorting_list.append([gene, pos, mutation]) #sort by position
	return [item[-1] for item in sorted(sorting_list)] #return original mutation names, but now sorted by position

def create_prevalences_input_table(mutations_of_interest, wdir, min_count, min_coverage, min_freq, mutation_count_file, mutation_coverage_file):
	'''
	import the PCA module which has genotype calling and
	filtering functions - outputs filtered AA tables, genotypes, and within sample allele frequencies to list only the samples
	that pass thresholds set ~3 cells above.
	'''

	subprocess.call(['mkdir', '-p', wdir])
	mutation_counts = get_mutation_counts(mutation_count_file)
	mutation_coverage = get_mutation_coverage(mutation_coverage_file)
	all_columns=mutation_counts.keys()
	mutation_dict={column[2]:column_number for column_number, column in enumerate(all_columns) if column[3]=='missense_variant'}

	#write only the desired mutations of interest to a new dataframe and output csv file
	mutations_of_interest=special_sort(mutations_of_interest)
	desired_column_numbers=[mutation_dict[mutation] for mutation in mutations_of_interest if mutation in mutation_dict]
	desired_columns=[all_columns[desired_column] for desired_column in desired_column_numbers]
	desired_mutation_counts=mutation_counts[desired_columns]
	desired_mutation_coverage=mutation_coverage[desired_columns]

	desired_mutation_counts.to_csv(os.path.join(wdir, 'desired_alternate_AA_table.csv'))
	desired_mutation_coverage.to_csv(os.path.join(wdir, 'desired_coverage_AA_table.csv'))

	gt_calls = PCA.call_genotypes(
		desired_mutation_counts, desired_mutation_coverage, min_count, min_coverage, min_freq
	)
	#print('options are', gt_calls.keys())

	# This step outputs samples that surpass the coverage threshold
	# as their original coverage value and resets samples not passing the threshold as '0'.
	# Note that even if there is insufficient coverage the resulting cell will be a '0' and not NaN
	filtered_coverage=gt_calls["filtered_mutation_coverage"]
	filtered_coverage.to_csv(os.path.join(wdir, "filtered_coverage_AA_table.csv"))

	# This step outputs samples containing alternate alleles that surpass the coverage and alternate thresholds
	# as their original values and samples not passing both thresholds as '0'. Note that even if there is insufficient
	# coverage to make a call the output for this table will still be listed as '0'
	filtered_mutation_counts = gt_calls["filtered_mutation_counts"]
	filtered_mutation_counts.to_csv(os.path.join(
			wdir, "filtered_alternate_AA_table.csv"))

	# This step outputs the within sample allele frequencies of each alternate allele. Unlike the tables above, this
	# one displays NaN if the coverage levels do not surpass the threshold (due to dividing alternate counts by 0)
	freq = gt_calls["wsaf"]
	freq.to_csv(os.path.join(
			wdir, "within_sample_allele_frequencies.csv"))

	# This step outputs a genotypes table, using the wsaf table as input. Cells can be NaN (coverage is below the
	# threshold), 0 (frequency of the mutation is less than or equal to min_freq), 1 (frequency is bigger than min_freq but
	# less than 100%) or 2 (mutation is found in all of the UMIs for the sample.
	genotypes = gt_calls["genotypes"]
	genotypes.to_csv(os.path.join(
			wdir, "filtered_genotypes_table.csv"))


	# The outputs of this table are the same as the genotypes table, except this table does not differentiate between
	# alternate alleles that are found in all UMIs of a sample and those that are found in some of the UMIs
	# (both '2' and '1' values from the genotypes table get reset to '1')
	prevalences_input = gt_calls["prevalences"]
	prevalences_input.to_csv(os.path.join(wdir, "prevalences_input_table.csv"))


	# To view the outputs of any given step, comment in the corresponding table below (or open the corresponding csv
	# file in a spreadsheet):
	# filtered_coverage.head()
	# filtered_mutation_counts.head()
	# genotypes.head()
	# prevalences_input.head()

def calculate_prevalences(wdir, metadata_file, sample_column, summary_column, mutations_of_interest):
	'''
	calculates final mutation prevalences for the mutations of interest and outputs a prevalence table. More specifically:
	prevalences_input_table = os.path.join(wdir, 'prevalences_input_table.csv')
	'''
	output_summary_table = os.path.join(wdir, 'prevalence_summary.tsv')
	prevalences_input_table = os.path.join(wdir, 'prevalences_input_table.csv')

	cap.calculate_prevalences(
		metadata_file,
		prevalences_input_table,
		mutations_of_interest,
		output_summary_table, 
		sample_column, 
		summary_column
	)

	prevalences=pd.read_csv(
		output_summary_table,
		header=list(range(1)),
		index_col=0, sep='\t'
	)

	return prevalences

def make_detail_graph(variant, wdir, zoom_level, latitude, longitude):
	df = pd.read_csv(os.path.join(wdir, "prevalence_summary.tsv"), sep='\t')
	summary_column = df.columns[0]
	if variant in list(df)[3:]:
		df["prevalence"] = df[variant].str.split(" ").str[0].astype(float)
		#df["prevalence_percent"] = df["prevalence"]*100+1
		#df["prevalence_log"]=np.log2(df["prevalence_percent"])
		max_prevalence = max(df["prevalence"].to_list())
		#max_prevalence=7
		#max_prevalence=1.0
		df["sample_size"] = (
			df[variant].str.split("/").str[1].str.replace(")", "").astype(float)
		)
		df = df[df["sample_size"] > 0]
		# df = df[df["Dataset"] == int(dataset)]

		fig = px.scatter_map(
			df,
			lat="Latitude",
			lon="Longitude",
			#color="prevalence_log",
			color="prevalence",
			size="sample_size",
			size_max=50,
			# color_continuous_scale='cividis',
			range_color=(0, max_prevalence), # use this to make the most prevalent mutation "pop" as bright yellow
			zoom=zoom_level,
			hover_name=summary_column,
			height=800,
			width=800,
			hover_data=["sample_size"],
			center={"lat": latitude, "lon": longitude},
			title=variant+'_'+summary_column+'_0-'+str(max_prevalence),
		)
	fig.update_layout(
		mapbox=dict(
			layers=[
				dict(
					sourcetype="geojson",
					source={"type": "FeatureCollection", "features": []},
					type="fill",
					color="blue",
					below="traces",
					visible=False  # Hide the layer by default
				)
			]
		)
	)
	return fig

def make_a_plot(sample_col, summary_col, variant, country, country_shortcuts, output_folder, metadata_table, mutations_of_interest, zoom, center_lat, center_long):
	calculate_prevalences(
		output_folder,
		metadata_table, 
		sample_col, 
		summary_col, 
		mutations_of_interest
	)
	if country != 'custom':
		zoom_level=country_shortcuts[country]['zoom']
		latitude=country_shortcuts[country]['lat']
		longitude=country_shortcuts[country]['lon']
	if country == 'custom':
		zoom_level = float(zoom)
		latitude = float(center_lat)
		longitude = float(center_long)
	fig = make_detail_graph(variant, output_folder, zoom_level, latitude, longitude)
	fig.show()
	
	title_string=f'{variant}_{summary_col}'
	fig.write_image(os.path.join(output_folder, title_string+'.svg'))
	fig.write_image(os.path.join(output_folder, title_string+'.png'))
	fig.write_html(os.path.join(output_folder, title_string+'.html'))

def get_countries_from_geojson():
	geojson_file = 'input/geojson_files/ne_adm0_10m.geojson'
	with open(geojson_file,'r') as f:
		adm0_json = json.load(f)
		adm0_set = set()
		for entry in adm0_json['features']:
			country = entry['properties']['NAME_EN']
			adm0_set.add(country)
		adm0_list = list(adm0_set)
		adm0_list.sort()
		# print(adm0_list)

	svg_country = widgets.Dropdown(
			options=adm0_list,
			description='Country:',
			disabled=False,
		)
	return svg_country

def create_static_maps(
	variant_of_interest,
	country_of_interest,
	graph_states_provinces,
	latitude_range,
	longitude_range,
	scale_factor,
	labels,
	annotation_font_size,
	title_text,
	title_size,
	wdir,
	summary_column):
	# prevalence parameters
	prevalence_df = pd.read_csv(os.path.join(wdir, "prevalence_summary.tsv"), sep='\t')
	prevalence_df["Prevalence"] = prevalence_df[variant_of_interest].str.split(" ").str[0].astype(float)
	prevalence_df["Population"] = prevalence_df[variant_of_interest].str.split("/").str[1].str.replace(")", "").astype(int)
	#prevalence_df['text'] = prevalence_df['Sites'] + ':' + prevalence_df['Population'].astype(str) + ', ' + prevalence_df['Prevalence'].astype(str)
	prevalence_df['text'] = prevalence_df[summary_column] + ':' + prevalence_df['Population'].astype(str) + ', ' + prevalence_df['Prevalence'].astype(str)
	prevalence_df['marker_size'] = prevalence_df['Population']*scale_factor
	prevalence_latitude = prevalence_df['Latitude']
	prevalence_longitude = prevalence_df['Longitude']

	fig = go.Figure()

	def plot_lakes():
		with open('input/geojson_files/ne_lakes_10m.geojson','r') as f:
			lakes_json = json.load(f)
			lakes_df = pd.json_normalize(lakes_json['features'])
			lakes_df['graphing_status'] = 1

		fig.add_choropleth(
			locations=lakes_df['properties.name'],
			z=lakes_df['graphing_status'],
			locationmode='geojson-id',
			geojson=lakes_json,
			featureidkey='properties.name',
			colorscale=['lightblue','lightblue'],
			marker_line_color='white',
			marker_line_width=0.1,
			showscale=False,
		)

	def plot_oceans():
		# don't need to do this if you just set plot background to blue
		with open('input/geojson_files/ne_ocean_10m.geojson','r') as f:
			ocean_json = json.load(f)
			ocean_df = pd.json_normalize(ocean_json['features'])
			ocean_df['graphing_status'] = 1

		fig.add_choropleth(
			locations=ocean_df['properties.featurecla'],
			z=ocean_df['graphing_status'],
			locationmode='geojson-id',
			geojson=ocean_json,
			featureidkey='properties.featurecla',
			colorscale=['lightblue','lightblue'],
			marker_line_color='white',
			marker_line_width=0.1,
			showscale=False,
		)

	def plot_adm0(adm0_country):
		with open('input/geojson_files/ne_adm0_10m.geojson','r') as f:
			adm0_json = json.load(f)
			adm0_df = pd.json_normalize(adm0_json['features'])
			adm0_df['graphing_status'] = 0
			adm0_df.loc[adm0_df['properties.NAME_EN'] == adm0_country, 'graphing_status'] = 1

		fig.add_choropleth(
			locations=adm0_df['properties.NAME_EN'],
			z=adm0_df['graphing_status'],
			locationmode='geojson-id',
			geojson=adm0_json,
			featureidkey='properties.NAME_EN',
			colorscale=['lightyellow','lightgray'],
			marker_line_color='black',
			marker_line_width=0.1,
			showscale=False,
		)

	def plot_adm1(adm1_country):
		with open('input/geojson_files/ne_adm1_10m.geojson', 'r') as f:
			adm1_json = json.load(f)
			adm1_json_filtered = adm1_json.copy()
			adm1_json_filtered['features'] = []
			for entry in adm1_json['features']:
				if entry['properties']['admin'] == adm1_country:
					adm1_json_filtered['features'].append(entry)
			adm1_df = pd.json_normalize(adm1_json_filtered['features'])
			adm1_df['graphing_status'] = 1
		if adm1_json_filtered['features'] != []:
			fig.add_choropleth(
				locations=adm1_df['properties.adm1_code'],
				z=adm1_df['graphing_status'],
				locationmode='geojson-id',
				geojson=adm1_json_filtered,
				featureidkey='properties.adm1_code',
				colorscale=['lightgray','lightgray'],
				marker_line_color='white',
				marker_line_width=0.5,
				showscale=False,
			)

	def plot_prevalence():
		fig.add_scattergeo(
			lon = prevalence_longitude,
			lat = prevalence_latitude,
			mode = 'markers+text',
			text = prevalence_df['text'],
			marker=dict(
				color = prevalence_df['Prevalence'],
				size = prevalence_df['marker_size'],
				sizemode = 'area',
				cmin = 0,
				cmax=1,
				showscale=True,
			),
			textposition='bottom right',
			textfont=dict(
				family="sans serif",
				size=annotation_font_size,
				color="black"
			),
			showlegend=False,
		)

	def add_labels():
		fig.add_scattergeo(
			lat=[label[1] for label in labels],
			lon=[label[2] for label in labels],
			text=[label[0] for label in labels],
			textposition="bottom center",
			textfont=dict(
			  family="sans serif",
			  size=annotation_font_size,
			  color="black"
			),
			mode='text',
			showlegend=False,
		)

	plot_adm0(country_of_interest)
	if graph_states_provinces:
		plot_adm1(country_of_interest)
	plot_lakes()
	plot_prevalence()
	add_labels()

	fig.update_geos(
		lataxis_range = latitude_range,
		lonaxis_range = longitude_range,
		visible=False,
		showrivers=True, rivercolor='lightblue',
		bgcolor='lightblue'
	)

	fig.update_layout(
		margin=dict(l=10, r=10, t=10, b=10),
		title=dict(text=title_text, font=dict(size=title_size), automargin=True, yref='paper')
	)

	subprocess.call(['mkdir', '-p', os.path.join(wdir, 'static_svg_files')])
	fig.write_image(os.path.join(wdir, 'static_svg_files', variant_of_interest+".svg"))
	fig.write_image(os.path.join(wdir, 'static_svg_files', variant_of_interest+".png"))
	# fig.write_image(os.path.join(wdir, 'static_svg_files', variant_of_interest+".html"))

	return fig