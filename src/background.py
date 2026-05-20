import sys
import calculate_prevalences as cap
import os
import pandas as pd
import plotly.express as px
import warnings
warnings.filterwarnings('ignore')
import ipywidgets as widgets
import yaml
import json
import plotly.graph_objs as go
import plotly.io as pio
import numpy as np
import PCA
import math
import subprocess
from functools import lru_cache
from urllib.request import urlopen

COUNTRIES_GEOJSON_URL = (
	"https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
	"geojson/ne_50m_admin_0_countries.geojson"
)
OCEAN_GEOJSON_URL = (
	"https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
	"geojson/ne_50m_ocean.geojson"
)
LOCAL_COUNTRIES_GEOJSON = 'input/geojson_files/ne_adm0_10m.geojson'

@lru_cache(maxsize=None)
def _load_geojson(geojson_file):
	with open(geojson_file, 'r') as f:
		return json.load(f)

@lru_cache(maxsize=1)
def _load_countries_geojson():
	try:
		with urlopen(COUNTRIES_GEOJSON_URL, timeout=10) as response:
			return json.load(response)
	except Exception:
		return _load_geojson(LOCAL_COUNTRIES_GEOJSON)

def _country_name(feature):
	properties = feature.get('properties', {})
	return (
		properties.get('NAME')
		or properties.get('NAME_EN')
		or properties.get('ADMIN')
		or properties.get('NAME_LONG')
	)

def _country_feature(country):
	for feature in _load_countries_geojson()['features']:
		if _country_name(feature) == country:
			return feature
	return None

@lru_cache(maxsize=1)
def _country_label_data():
	labels = []
	for feature in _load_countries_geojson()['features']:
		country = _country_name(feature)
		if not country:
			continue
		properties = feature.get('properties', {})
		longitude = properties.get('LABEL_X')
		latitude = properties.get('LABEL_Y')
		if longitude is None or latitude is None:
			latitudes, longitudes = [], []
			_collect_geojson_coordinates(feature['geometry'], latitudes, longitudes)
			if not latitudes or not longitudes:
				continue
			longitude = (min(longitudes) + max(longitudes)) / 2
			latitude = (min(latitudes) + max(latitudes)) / 2
		labels.append((country, latitude, longitude))
	return labels

def _collect_geojson_coordinates(geometry, latitudes, longitudes):
	geometry_type = geometry.get('type')
	coordinates = geometry.get('coordinates', [])
	if geometry_type == 'Point':
		longitudes.append(coordinates[0])
		latitudes.append(coordinates[1])
	elif geometry_type in {'MultiPoint', 'LineString'}:
		for longitude, latitude in coordinates:
			longitudes.append(longitude)
			latitudes.append(latitude)
	elif geometry_type in {'MultiLineString', 'Polygon'}:
		for line in coordinates:
			for longitude, latitude in line:
				longitudes.append(longitude)
				latitudes.append(latitude)
	elif geometry_type == 'MultiPolygon':
		for polygon in coordinates:
			for line in polygon:
				for longitude, latitude in line:
					longitudes.append(longitude)
					latitudes.append(latitude)
	elif geometry_type == 'GeometryCollection':
		for child_geometry in geometry.get('geometries', []):
			_collect_geojson_coordinates(child_geometry, latitudes, longitudes)

def _feature_intersects_bounds(feature, latitude_range, longitude_range):
	latitudes, longitudes = [], []
	_collect_geojson_coordinates(feature['geometry'], latitudes, longitudes)
	if not latitudes or not longitudes:
		return False
	return (
		min(latitudes) <= max(latitude_range)
		and max(latitudes) >= min(latitude_range)
		and min(longitudes) <= max(longitude_range)
		and max(longitudes) >= min(longitude_range)
	)

@lru_cache(maxsize=None)
def _load_geojson_in_bounds(geojson_file, latitude_range_key, longitude_range_key):
	geojson = _load_geojson(geojson_file)
	return _filter_geojson_in_bounds(geojson, latitude_range_key, longitude_range_key)

def _filter_geojson_in_bounds(geojson, latitude_range_key, longitude_range_key):
	latitude_range = list(latitude_range_key)
	longitude_range = list(longitude_range_key)
	return {
		**geojson,
		'features': [
			feature for feature in geojson['features']
			if _feature_intersects_bounds(feature, latitude_range, longitude_range)
		]
	}

def _range_cache_key(axis_range):
	return tuple(round(float(value), 6) for value in axis_range)

def _padded_bounds_from_center(latitude, longitude, zoom_level, padding_factor=1.5):
	span = max(0.5, 2 ** (8.5 - zoom_level))
	latitude_range = [
		max(latitude - span * padding_factor / 2, -90),
		min(latitude + span * padding_factor / 2, 90),
	]
	longitude_range = [
		max(longitude - span * padding_factor / 2, -180),
		min(longitude + span * padding_factor / 2, 180),
	]
	return _range_cache_key(latitude_range), _range_cache_key(longitude_range)

def get_country_geo_ranges(country, padding_fraction=0.08, minimum_padding=0.25):
	adm0_json = _load_countries_geojson()
	for entry in adm0_json['features']:
		if _country_name(entry) == country:
			latitudes, longitudes = [], []
			_collect_geojson_coordinates(entry['geometry'], latitudes, longitudes)
			if not latitudes or not longitudes:
				break
			latitude_padding = max((max(latitudes) - min(latitudes)) * padding_fraction, minimum_padding)
			longitude_padding = max((max(longitudes) - min(longitudes)) * padding_fraction, minimum_padding)
			return (
				[max(min(latitudes) - latitude_padding, -90), min(max(latitudes) + latitude_padding, 90)],
				[max(min(longitudes) - longitude_padding, -180), min(max(longitudes) + longitude_padding, 180)]
			)
	raise ValueError(f'Could not find country {country!r} in countries GeoJSON')

def _zoom_for_bounds(latitude_range, longitude_range):
	latitude_span = max(latitude_range) - min(latitude_range)
	longitude_span = max(longitude_range) - min(longitude_range)
	largest_span = max(latitude_span, longitude_span, 0.01)
	return max(1, min(9, 8.5 - math.log2(largest_span)))

def get_country_map_view(country):
	latitude_range, longitude_range = get_country_geo_ranges(country)
	return {
		'lat': sum(latitude_range) / 2,
		'lon': sum(longitude_range) / 2,
		'zoom': _zoom_for_bounds(latitude_range, longitude_range),
		'bounds': {
			'south': min(latitude_range),
			'north': max(latitude_range),
			'west': min(longitude_range),
			'east': max(longitude_range),
		},
	}

def _should_auto_range(axis_range):
	return axis_range is None or axis_range == 'auto'

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

def get_mutation_counts(mutation_count_file, mutation_coverage_file):
	mutation_counts = pd.read_csv(mutation_count_file,
								  header=list(range(6)),
								  index_col=0)
	mutation_coverage = pd.read_csv(mutation_coverage_file,
									index_col=0,
									header=list(range(6)))
	all_columns=list(mutation_counts.keys())
	mutation_dict={column[2]:column_number for column_number, column in enumerate(all_columns) if column[3]=='missense_variant'}
	filtered_mutations=special_sort(list(mutation_dict.keys()))
	return all_columns, mutation_dict, filtered_mutations, mutation_counts, mutation_coverage

def get_metadata_columns(metadata_table, separator='\t'):
	separator = cap.detect_separator(metadata_table, separator)
	metadata = pd.read_csv(metadata_table, sep=separator)
	columns = list(metadata.keys())
	return columns

def generate_country_dropdown():
	countries_json = _load_countries_geojson()
	country_names = sorted(
		country for country in {_country_name(feature) for feature in countries_json['features']}
		if country
	)
	country_shortcuts = {country: get_country_map_view(country) for country in country_names}

	country = widgets.Dropdown(
		options=country_names,
		description='Country:',
		disabled=False,
	)
	if 'Uganda' in country_names:
		country.value = 'Uganda'
	return country, country_shortcuts

# RUN

def create_prevalences_input_table(mutations_of_interest, mutation_dict, all_columns, mutation_counts, mutation_coverage, wdir, min_count, min_coverage, min_freq):
	'''
	import the PCA module which has genotype calling and
	filtering functions - outputs filtered AA tables, genotypes, and within sample allele frequencies to list only the samples
	that pass thresholds set ~3 cells above.
	'''

	subprocess.call(['mkdir', '-p', wdir])

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
	freq.head()

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

def calculate_prevalences(wdir, metadata_file, sample_counts_file, sample_column, summary_column, mutations_of_interest, separator='\t'):
	'''
	calculates final mutation prevalences for the mutations of interest and outputs a prevalence table. More specifically:
	prevalences_input_table = os.path.join(wdir, 'prevalences_input_table.csv')
	'''
	output_summary_table = os.path.join(wdir, 'prevalence_summary.tsv')
	prevalences_input_table = os.path.join(wdir, sample_counts_file)
	# print('sample column is', sample_column)
	# print('summary column is', summary_column)
	first_line=open(metadata_file).readline()
	# print('first line of metadata file is', first_line)

	cap.calculate_prevalences(metadata_file,
						  prevalences_input_table,
						  mutations_of_interest,
						  output_summary_table, sample_column, summary_column, separator)
	prevalences=pd.read_csv(output_summary_table,
							  header=list(range(1)),
							  index_col=0, sep='\t')
	return prevalences


def add_countries_only_layer(
	fig,
	chosen_country=None,
	overlay_water=False,
	highlight_chosen_country=True,
	country_bounds=None,
	latitude=None,
	longitude=None,
	zoom_level=None):
	map_layers = [
		{
			"sourcetype": "geojson",
			"source": {
				"type": "FeatureCollection",
				"features": [
					{
						"type": "Feature",
						"geometry": {
							"type": "Polygon",
							"coordinates": [[
								[-180, -90],
								[180, -90],
								[180, 90],
								[-180, 90],
								[-180, -90],
							]],
						},
						"properties": {},
					}
				],
			},
			"type": "fill",
			"color": "lightblue",
			"below": "traces",
		}
	]
	if overlay_water and country_bounds is not None:
		latitude_span = country_bounds['north'] - country_bounds['south']
		longitude_span = country_bounds['east'] - country_bounds['west']
		latitude_padding = max(latitude_span * 0.15, 0.25)
		longitude_padding = max(longitude_span * 0.15, 0.25)
		latitude_range_key = _range_cache_key([
			max(country_bounds['south'] - latitude_padding, -90),
			min(country_bounds['north'] + latitude_padding, 90),
		])
		longitude_range_key = _range_cache_key([
			max(country_bounds['west'] - longitude_padding, -180),
			min(country_bounds['east'] + longitude_padding, 180),
		])
	elif overlay_water and latitude is not None and longitude is not None and zoom_level is not None:
		latitude_range_key, longitude_range_key = _padded_bounds_from_center(latitude, longitude, zoom_level, padding_factor=0.8)
	else:
		latitude_range_key = None
		longitude_range_key = None
	if overlay_water and latitude_range_key is not None and longitude_range_key is not None:
		water_layers = []
		lake_layers = [
			{
				"sourcetype": "geojson",
				"source": _load_geojson_in_bounds(
					'input/geojson_files/ne_lakes_10m.geojson',
					latitude_range_key,
					longitude_range_key
				),
				"type": "fill",
				"color": "lightblue",
				"below": "traces",
			},
		]
	else:
		water_layers = []
		lake_layers = []
	map_layers.extend(water_layers)
	map_layers.append(
		{
			"sourcetype": "geojson",
			"source": COUNTRIES_GEOJSON_URL,
			"type": "fill",
			"color": "lightgray",
			"below": "traces",
		}
	)
	chosen_feature = _country_feature(chosen_country) if highlight_chosen_country and chosen_country is not None else None
	if chosen_feature is not None:
		map_layers.append(
			{
				"sourcetype": "geojson",
				"source": {"type": "FeatureCollection", "features": [chosen_feature]},
				"type": "fill",
				"color": "white",
				"below": "traces",
			}
		)
	map_layers.extend(lake_layers)
	map_layers.append(
		{
			"sourcetype": "geojson",
			"source": COUNTRIES_GEOJSON_URL,
			"type": "line",
			"color": "black",
			"line": {"width": 1},
			"below": "traces",
		}
	)
	fig.update_layout(
		map_style="white-bg",
		map_layers=map_layers
	)
	country_labels = _country_label_data()
	fig.add_trace(
		go.Scattermap(
			lon=[label[2] for label in country_labels],
			lat=[label[1] for label in country_labels],
			text=[label[0] for label in country_labels],
			mode="text",
			textfont=dict(size=8, color="dimgray"),
			hoverinfo="skip",
			showlegend=False,
		)
	)
	fig.data = (fig.data[-1],) + fig.data[:-1]
	return fig

def _detail_prevalence_df(variant, summary_column, wdir):
	df = pd.read_csv(os.path.join(wdir, "prevalence_summary.tsv"), sep='\t')
	if variant not in list(df)[3:]:
		raise ValueError(f'{variant!r} is not found in prevalence_summary.tsv')
	df["prevalence"] = df[variant].str.split(" ").str[0].astype(float)
	df["sample_size"] = (
		df[variant].str.split("/").str[1].str.replace(")", "").astype(float)
	)
	df = df[df["sample_size"] > 0].copy()
	df["point_label"] = (
		df[summary_column].astype(str)
		+ ", "
		+ df["sample_size"].astype("Int64").astype(str)
		+ ", "
		+ df["prevalence"].astype(str)
	)
	return df

def make_detail_graph(
	variant,
	summary_column,
	wdir,
	zoom_level,
	latitude,
	longitude,
	country_bounds=None,
	chosen_country=None,
	overlay_water=False,
	highlight_chosen_country=True,
	use_geojson=True,
	use_country_geojson=None):
	if use_country_geojson is not None:
		use_geojson = use_country_geojson
	df = _detail_prevalence_df(variant, summary_column, wdir)
		#df["prevalence_percent"] = df["prevalence"]*100+1
		#df["prevalence_log"]=np.log2(df["prevalence_percent"])
	max_prevalence = max(df["prevalence"].to_list())
		#max_prevalence=7
		#max_prevalence=1.0
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
		text="point_label",
		height=800,
		width=800,
		hover_data=["sample_size"],
		center={"lat": latitude, "lon": longitude},
		title=variant+'_'+summary_column+'_0-'+str(max_prevalence),
	)
	fig.update_traces(
		mode="markers+text",
		textposition="middle right",
		textfont=dict(size=10, color="black"),
	)
	if use_geojson:
		fig = add_countries_only_layer(
				fig,
				chosen_country,
				overlay_water=overlay_water,
				highlight_chosen_country=highlight_chosen_country,
				country_bounds=country_bounds,
			latitude=latitude,
			longitude=longitude,
			zoom_level=zoom_level
		)
	else:
		fig.update_layout(map_style="open-street-map")
	return fig

def make_detail_graph_vector(
	variant,
	summary_column,
	wdir,
	country_bounds=None,
	chosen_country=None,
	overlay_water=False,
	highlight_chosen_country=True,
	width=800,
	height=800):
	df = _detail_prevalence_df(variant, summary_column, wdir)
	max_prevalence = max(df["prevalence"].to_list())
	countries_json = _load_countries_geojson()
	country_names = [_country_name(feature) for feature in countries_json['features']]

	fig = go.Figure()
	fig.add_choropleth(
		locations=country_names,
		z=[0] * len(country_names),
		locationmode='geojson-id',
		geojson=countries_json,
		featureidkey='properties.NAME',
		colorscale=['lightgray', 'lightgray'],
		marker_line_color='black',
		marker_line_width=0.5,
		showscale=False,
	)

	chosen_feature = _country_feature(chosen_country) if highlight_chosen_country and chosen_country is not None else None
	if chosen_feature is not None:
		fig.add_choropleth(
			locations=[_country_name(chosen_feature)],
			z=[1],
			locationmode='geojson-id',
			geojson={"type": "FeatureCollection", "features": [chosen_feature]},
			featureidkey='properties.NAME',
			colorscale=['white', 'white'],
			marker_line_color='black',
			marker_line_width=0.5,
			showscale=False,
		)

	if overlay_water and country_bounds is not None:
		latitude_span = country_bounds['north'] - country_bounds['south']
		longitude_span = country_bounds['east'] - country_bounds['west']
		latitude_padding = max(latitude_span * 0.15, 0.25)
		longitude_padding = max(longitude_span * 0.15, 0.25)
		latitude_range_key = _range_cache_key([
			max(country_bounds['south'] - latitude_padding, -90),
			min(country_bounds['north'] + latitude_padding, 90),
		])
		longitude_range_key = _range_cache_key([
			max(country_bounds['west'] - longitude_padding, -180),
			min(country_bounds['east'] + longitude_padding, 180),
		])
		lakes_json = _load_geojson_in_bounds(
			'input/geojson_files/ne_lakes_10m.geojson',
			latitude_range_key,
			longitude_range_key
		)
		if lakes_json['features']:
			lake_names = [
				feature.get('properties', {}).get('name') or str(index)
				for index, feature in enumerate(lakes_json['features'])
			]
			for lake_name, feature in zip(lake_names, lakes_json['features']):
				feature.setdefault('properties', {})['plot_id'] = lake_name
			fig.add_choropleth(
				locations=lake_names,
				z=[0] * len(lake_names),
				locationmode='geojson-id',
				geojson=lakes_json,
				featureidkey='properties.plot_id',
				colorscale=['lightblue', 'lightblue'],
				marker_line_color='lightblue',
				marker_line_width=0,
				showscale=False,
			)

	country_labels = _country_label_data()
	fig.add_scattergeo(
		lon=[label[2] for label in country_labels],
		lat=[label[1] for label in country_labels],
		text=[label[0] for label in country_labels],
		mode='text',
		textfont=dict(size=8, color='dimgray'),
		hoverinfo='skip',
		showlegend=False,
	)

	fig.add_scattergeo(
		lon=df["Longitude"],
		lat=df["Latitude"],
		text=df["point_label"],
		mode='markers+text',
		textposition='middle right',
		textfont=dict(size=10, color='black'),
		marker=dict(
			color=df["prevalence"],
			size=df["sample_size"],
			sizemode='area',
			cmin=0,
			cmax=max_prevalence,
			showscale=True,
			colorbar=dict(title='prevalence'),
		),
		hovertext=df[summary_column],
		showlegend=False,
	)

	if country_bounds is not None:
		lat_range = [country_bounds['south'], country_bounds['north']]
		lon_range = [country_bounds['west'], country_bounds['east']]
	else:
		lat_range = [df["Latitude"].min() - 1, df["Latitude"].max() + 1]
		lon_range = [df["Longitude"].min() - 1, df["Longitude"].max() + 1]
	fig.update_geos(
		lataxis_range=lat_range,
		lonaxis_range=lon_range,
		visible=False,
		bgcolor='lightblue',
	)
	fig.update_layout(
		width=width,
		height=height,
		margin=dict(l=10, r=10, t=40, b=10),
		title=variant+'_'+summary_column+'_0-'+str(max_prevalence),
	)
	return fig

def write_detail_graph_files(
	fig,
	output_folder,
	title_string,
	variant,
	summary_column,
	country_bounds=None,
	chosen_country=None,
	overlay_water=False,
	highlight_chosen_country=True):
	write_detail_graph_vector_files(
		variant,
		summary_column,
		output_folder,
		title_string,
		country_bounds=country_bounds,
		chosen_country=chosen_country,
		overlay_water=overlay_water,
		highlight_chosen_country=highlight_chosen_country,
		width=fig.layout.width or 800,
		height=fig.layout.height or 800,
	)
	fig.write_html(os.path.join(output_folder, title_string+'.html'))

def write_detail_graph_vector_files(
	variant,
	summary_column,
	wdir,
	title_string,
	country_bounds=None,
	chosen_country=None,
	overlay_water=False,
	highlight_chosen_country=True,
	width=800,
	height=800):
	import geopandas as gpd
	import matplotlib.pyplot as plt

	df = _detail_prevalence_df(variant, summary_column, wdir)
	countries = gpd.read_file(LOCAL_COUNTRIES_GEOJSON)
	if country_bounds is not None:
		west = country_bounds['west']
		east = country_bounds['east']
		south = country_bounds['south']
		north = country_bounds['north']
	else:
		west = df["Longitude"].min() - 1
		east = df["Longitude"].max() + 1
		south = df["Latitude"].min() - 1
		north = df["Latitude"].max() + 1
	countries = countries.cx[west:east, south:north]

	fig_size = (width / 100, height / 100)
	mpl_fig, ax = plt.subplots(figsize=fig_size)
	ax.set_facecolor('lightblue')
	countries.plot(ax=ax, color='lightgray', edgecolor='none', linewidth=0)
	if highlight_chosen_country and chosen_country is not None:
		name_columns = [column for column in ['NAME', 'NAME_EN', 'ADMIN', 'NAME_LONG'] if column in countries.columns]
		if name_columns:
			mask = False
			for column in name_columns:
				mask = mask | (countries[column] == chosen_country)
			countries[mask].plot(ax=ax, color='white', edgecolor='none', linewidth=0)
	if overlay_water:
		lakes = gpd.read_file('input/geojson_files/ne_lakes_10m.geojson')
		lakes = lakes.cx[west:east, south:north]
		if not lakes.empty:
			lakes.plot(ax=ax, color='lightblue', edgecolor='lightblue', linewidth=0)
	countries.boundary.plot(ax=ax, color='black', linewidth=0.5)

	for country_name, label_latitude, label_longitude in _country_label_data():
		if west <= label_longitude <= east and south <= label_latitude <= north:
			ax.text(label_longitude, label_latitude, country_name, fontsize=8, color='dimgray', ha='center', va='center')

	sizes = df["sample_size"].clip(lower=1) * 8
	points = ax.scatter(
		df["Longitude"],
		df["Latitude"],
		c=df["prevalence"],
		s=sizes,
		cmap='viridis',
		vmin=0,
		vmax=max(df["prevalence"].max(), 1e-9),
		edgecolors='black',
		linewidths=0.5,
	)
	for _, row in df.iterrows():
		ax.annotate(
			row["point_label"],
			(row["Longitude"], row["Latitude"]),
			textcoords='offset points',
			xytext=(8, 0),
			ha='left',
			va='center',
			fontsize=10,
			color='black',
		)
	mpl_fig.colorbar(points, ax=ax, label='prevalence', fraction=0.035, pad=0.02)
	ax.set_xlim(west, east)
	ax.set_ylim(south, north)
	ax.set_axis_off()
	ax.set_title(variant+'_'+summary_column+'_0-'+str(max(df["prevalence"].to_list())))
	mpl_fig.tight_layout()
	mpl_fig.savefig(os.path.join(wdir, title_string+'.svg'), format='svg')
	mpl_fig.savefig(os.path.join(wdir, title_string+'.png'), dpi=300)
	plt.close(mpl_fig)

def get_countries_from_geojson():
	geojson_file = 'input/geojson_files/ne_adm0_10m.geojson'
	adm0_json = _load_geojson(geojson_file)
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
	summary_column,
	write_files=False):
	if _should_auto_range(latitude_range) or _should_auto_range(longitude_range):
		auto_latitude_range, auto_longitude_range = get_country_geo_ranges(country_of_interest)
		if _should_auto_range(latitude_range):
			latitude_range = auto_latitude_range
		if _should_auto_range(longitude_range):
			longitude_range = auto_longitude_range
	latitude_range_key = _range_cache_key(latitude_range)
	longitude_range_key = _range_cache_key(longitude_range)

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
		lakes_json = _load_geojson_in_bounds(
			'input/geojson_files/ne_lakes_10m.geojson',
			latitude_range_key,
			longitude_range_key
		)
		if not lakes_json['features']:
			return
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
		ocean_json = _load_geojson('input/geojson_files/ne_ocean_10m.geojson')
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
		adm0_json = _load_geojson_in_bounds(
			'input/geojson_files/ne_adm0_10m.geojson',
			latitude_range_key,
			longitude_range_key
		)
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
		adm1_json = _load_geojson('input/geojson_files/ne_adm1_10m.geojson')
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
		width=800, height=600,
		margin=dict(l=10, r=10, t=10, b=10),
		title=dict(text=title_text, font=dict(size=title_size), automargin=True, yref='paper')
	)

	if write_files:
		static_map_folder = os.path.join(wdir, 'static_svg_files')
		os.makedirs(static_map_folder, exist_ok=True)
		fig.write_image(os.path.join(static_map_folder, variant_of_interest+".svg"))
		fig.write_image(os.path.join(static_map_folder, variant_of_interest+".png"))
		fig.write_html(os.path.join(static_map_folder, variant_of_interest+".html"))

	return fig
