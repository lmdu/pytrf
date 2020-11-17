#include "tandem.h"

PyObject *find_short_tandem_repeats(PyObject *self, PyObject *args, PyObject *kwargs) {
	char *seq;
	Py_ssize_t size;

	//minimum repeat array
	uint32_t minreps[7] = {0, 0, 0, 0, 0, 0, 0};

	// current slide position
	uint64_t i;

	//motif length
	uint8_t j;

	// tandem repeat start location
	uint64_t start;

	// tandem repeat length
	uint32_t replen;

	//tandem repeat repeats
	uint32_t repeats;

	//motif sequence
	char motif[7];

	char *keywords[] = {"seq", "minreps", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#(iiiiii)", keywords, &seq, &size, &minreps[1],
		&minreps[2], &minreps[3], &minreps[4], &minreps[5], &minreps[6])) {
		return NULL;
	}

	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;

	for (i = 0; i < size; ++i) {
		if (seq[i] == 78) {
			continue;
		}

		for (j = 1; j <= 6; ++j) {
			start = i;

			while (seq[i] == seq[i+j]) {
				++i;
			}

			replen = i + j - start;
			repeats = replen/j;

			if (repeats >= minreps[j]) {
				memcpy(motif, seq+start, j);
				motif[j] = '\0';

				i += j;

				if (motif[0] != 78) {
					tmp = Py_BuildValue("(siiiii)", motif, j, repeats, start+1, i, replen);
					PyList_Append(ssrs, tmp);
					Py_DECREF(tmp);
				}

				j = 0;
			} else {
				i = start;
			}
		}
	}
	return ssrs;
}

