from zipfile import ZipFile
import xml.etree.ElementTree as et
from io import BytesIO
import numpy as np


class UNICORNZipExport():
	'''
	Read UNICORN 7 bulk export (ZIP) files

	Exposes the data in the UNICORN file at the path that the class was
	constructed with. All data is stored as parallel arrays.

	data_order
		Whether the volume measurement or read-out comes first in the
		data points, e.g.:

		(0.2 mL, 0.35 mAU)

		or

		(0.35 mAU, 0.2 mL).

		The default is volume first, which can be changed by providing
		the class constructor with a reversed `tuple`.

	ignore_missing_data
		Set this to `True` if the class should just keep going if
		expected data is not present.

	fractions
		The volume at which each fraction begins and corresponding code.

	injections
		The volume and name of injections of material.

	messages
		The volume and description of system events that happened
		during the run.

	readout
		Data concerning the different readouts from the sensors of the
		ÄKTA machine.
	'''

	data_order = ('volume', 'amplitude')
	fractions = []
	injections = []
	messages = []
	readout = {}

	def __init__(
		self,
		file_name,
		data_order=('volume', 'amplitude'),
		ignore_missing_data=False
	):

		if 'volume' in data_order and 'amplitude' in data_order:
			self._do = data_order
		else:
			raise ValueError('invalid data order')

		self._load_data(file_name, ignore_missing_data)

	def _load_data(self, f_name, ignore_missing_data):
		with ZipFile(f_name) as ar:
			metadata = 'Chrom.1.Xml'
			with ar.open(metadata) as mtd:
				root = et.parse(mtd)

				events = self._get_events(root)
				try:
					self.fractions = events['Fraction']
					self.injections = events['Injection']
					self.messages = events['Logbook']
				except KeyError as e:
					if ignore_missing_data:
						pass
					else:
						raise ValueError(
								f'exported file "{f_name}" does not include "{e.args[0]}" events'
							)

				self.readout = self._get_ro(ar, root)

	class _ReadOut():
		'''
		Data for a single ÄKTA read-out

		measurements
			Volume and amplitude data for this read-out.

		units
			The units of the above data.

		wavelength
			For optical read-outs, the wavelength of the light used
			in nanometers. `None` if not applicable.
		'''

		measurements = []
		units = ()
		wavelength = None

	def _get_ro(self, archive, xml_root):
		ros = {}
		for curve in xml_root.iter('Curve'):
			ro = self._ReadOut()

			[f_name] = curve.iter('BinaryCurvePointsFileName')
			ro.measurements = self._decode_ro_buf(archive.read(f_name.text))

			ro.units = self._get_ro_units(curve)

			cn = curve.find('Name')
			name, wl = self._correct_curve_name(cn.text)
			ro.wavelength = wl

			ros[name] = ro
		return ros

	def _get_ro_units(self, ro_xml):
		unit_full_names = (n.capitalize() + 'Unit' for n in self._do)
		return (ro_xml.find(fn).text for fn in unit_full_names)

	@staticmethod
	def _correct_curve_name(n):
		if n.startswith('UV') and '_' in n:
			nice, wl = n.split('_')
			return nice, int(wl)
		else:
			better_names = {
			'Cond'			: 'Conductance',
			'% Cond'		: 'Conductance (%)',
			'Conc B'		: 'Concentration B',
			'PreC pressure'	: 'Pre-column pressure',
			'Cond temp'		: 'Conductance temperature',	# ???
			}
			return better_names.get(n, n), None

	def _get_events(self, xml_root):
		event_per_type = {}
		for ec in xml_root.iter('EventCurve'):
			events = [
				(
					float(e.find('EventVolume').text),
					e.find('EventText').text,
				)
				for e in ec.iter('Event')
			]
			nice = tuple(list(r) for r in zip(*events))
			if self._do[0] != 'volume':
				nice = reverse(nice)

			type = ec.get('EventCurveType')

			event_per_type[type] = nice
		return event_per_type

	def _decode_ro_buf(self, buf):
		stop = self._get_padding_start(buf)
		# `ZipFile` wants to `seek()`, but does not allow specifing a
		# byte offset, so camouflage the bytes object as a buffer
		with ZipFile(BytesIO(buf[:stop])) as ar:
			dim_f_names = [
				f'CoordinateData.{d.capitalize()}s'
				for d in self._do
			]

			dim_data = [
				self._decode_binary_data(ar.read(n))
				for n in dim_f_names 
			]

		return dim_data

	@staticmethod
	def _get_padding_start(buf):
		eocd_start = buf.rindex(b'PK\x05\x06')
		eocd_len = 22	# assuming no ZIP comment
		zip_end = eocd_start + eocd_len

		return zip_end

	@staticmethod
	def _decode_binary_data(buf):
		prefix_len = 3 + 6 * 4	# not clear why
		float_len = 4
		buf_end = prefix_len % float_len - float_len	# to align buffer
		floats = np.frombuffer(
			buf[:buf_end],
			dtype=np.dtype(f'<f{float_len}'),
			offset=prefix_len,
		)

		return floats
