/*
* This is utility class to transpose FileVector files in binary format, so there will not be need in
advanced text tools behaviour.
*/

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "frutil.h"
#include "FileVector.h"
#include "Transposer.h"

using namespace std;

void Transposer::process(string filename) {
    process(filename, string(""), false );
}


void Transposer::process(string filename, string destFileName,
                         bool forceOverwrite) {
    FileVector* src_fv = new FileVector(filename,1);
    unsigned long src_nvars = src_fv->getNumVariables();
    unsigned long src_nobss = src_fv->getNumObservations();
    unsigned int data_size = src_fv->getElementSize();

    string dest_file_name;
    string src_data_file_name;
    string dest_data_file_name;

    if (destFileName == "") {
        // legacy
        dest_file_name = extract_base_file_name(filename) + "_transposed";
        src_data_file_name = extract_base_file_name(filename) +
            FILEVECTOR_DATA_FILE_SUFFIX;
        dest_data_file_name = extract_base_file_name(filename) + "_transposed" +
            FILEVECTOR_DATA_FILE_SUFFIX;
    } else {
        dest_file_name = destFileName;
        src_data_file_name = filename + FILEVECTOR_DATA_FILE_SUFFIX;
        dest_data_file_name = destFileName + FILEVECTOR_DATA_FILE_SUFFIX;
    }


    if (!forceOverwrite && headerOrDataExists(dest_file_name)) {
        errorLog << "File already exists: " << dest_file_name
                 << endl << errorExit;
    }

    initializeEmptyFile(dest_file_name, src_fv->getNumObservations(),
                        src_fv->getNumVariables(), src_fv->getElementType(),
                        true);

    FileVector* dest_fv = new FileVector(dest_file_name,1);
    dbg << "Copying var/obs names...";
    write_var_obs_names(src_fv, dest_fv);

    delete src_fv;
    delete dest_fv;
    dbg << "done" << endl;

    copy_data(src_data_file_name, dest_data_file_name, src_nvars, src_nobss,
              data_size);
    dbg << "done"<< endl;
}


void Transposer::write_var_obs_names(FileVector *src_fv, FileVector *dest_fv) {
    // copy observations and variables names
    for (unsigned long i = 0 ; i < src_fv->getNumVariables(); i++ )
        dest_fv->writeObservationName( i, src_fv->readVariableName(i) );

    for (unsigned long i = 0 ; i < src_fv->getNumObservations(); i++ )
        dest_fv->writeVariableName( i, src_fv->readObservationName(i) );
}


void Transposer::copy_data(string src_data_file_name, string dest_data_file_name,
                           unsigned long src_nvars, unsigned long src_nobss,
                           unsigned int data_size) {
    dbg << "Copying data..." << src_nobss << "x" << src_nvars << endl;

    unsigned long obs_pages = src_nobss / square_size;
    if (src_nobss % square_size > 0) obs_pages++;

    unsigned long var_pages = src_nvars / square_size;
    if (src_nvars % square_size > 0) var_pages++;

    ifstream * src_stream = new ifstream();
    src_stream->open(src_data_file_name.c_str(),ios::in | ios::binary);

    ofstream * dest_stream = new ofstream;
    dest_stream->open(dest_data_file_name.c_str(),ios::out | ios::binary);

    for (unsigned long i = 0; i < var_pages; i++) {
        for (unsigned long j = 0; j< obs_pages; j++) {
            unsigned long obs_length = square_size;
            if ((j + 1 ) * square_size > src_nobss)
                obs_length = src_nobss % square_size;

            unsigned long var_length = square_size;
            if ((i + 1 ) * square_size > src_nvars)
                var_length = src_nvars % square_size;

            char * data_part =
                new (nothrow) char[var_length*obs_length*data_size];
            if (!data_part)
                errorLog << "can not allocate memory for data_part"
                         << errorExit;

            char * data_part_transposed =
                new (nothrow) char[var_length*obs_length*data_size];
            if (!data_part_transposed)
                errorLog << "can not allocate memory for data_part_transposed"
                         << errorExit;

            read_part(src_stream, data_part, j * square_size ,
                      obs_length, i * square_size , var_length,
                      data_size, src_nobss);
            transpose_part(data_part, data_part_transposed, obs_length,
                           var_length, data_size);
            write_part(dest_stream, data_part_transposed, i * square_size,
                       var_length, j * square_size, obs_length,
                       data_size, src_nvars);

            delete[] data_part;
            delete[] data_part_transposed;
        }
        dbg << endl;
    }

    src_stream->close();
    delete src_stream;
    dest_stream->close();
    delete dest_stream;

    dbg << "data written" << endl;
}


/*
 * read next piece of data with size = obs_length x var_length,
 * starting from var_start, obs_start coordinates
 */
void Transposer::read_part(ifstream * src_stream, char * data_part,
                           unsigned long obs_start, unsigned long obs_length,
                           unsigned long var_start, unsigned long var_length,
                           unsigned int  data_size,
                           unsigned long src_obs_length) {
    for (unsigned long i = 0; i < var_length; i++) {
        //seek to the beginning of the next var
        unsigned long read_pos = (var_start + i ) * src_obs_length + obs_start;
        src_stream->seekg(read_pos * data_size, ios::beg);
        //read next var to input buffer
        src_stream->read(data_part + ( i * obs_length * data_size ),
                         obs_length * data_size );
    }
}


/*
 * write next piece of transposed data with size = obs_length' x var_length'
 */
void Transposer::write_part(ofstream * dest_stream, char * data_part_transposed,
                            unsigned long obs_start, unsigned long obs_length,
                            unsigned long var_start, unsigned long var_length,
                            unsigned int  data_size,
                            unsigned long dest_obs_length) {
    for (unsigned long i = 0; i < var_length; i++) {
        // seek to the beginning of the next var
        unsigned long write_pos = (var_start + i ) * dest_obs_length + obs_start;
        dest_stream->seekp(write_pos * data_size, ios::beg);
        // write next piece of var to file
        dest_stream->write(data_part_transposed + (i * obs_length * data_size),
                           obs_length * data_size );
    }
}


/*
 * transpose piece of data to write to the new file.
 * original axb matrix flipped to bxa matrix.
 */
void Transposer::transpose_part(void * data_part, void * data_part_transposed,
                                unsigned long obs_length,
                                unsigned long var_length,
                                unsigned int data_size) {
    for (unsigned long i = 0; i < var_length; i++) {
        for(unsigned long j = 0; j < obs_length; j++) {
            int from_pos = (i * obs_length + j ) * data_size;
            int to_pos = ( j * var_length  + i ) * data_size;
            memcpy((char*)data_part_transposed + to_pos,
                   (char*)data_part + from_pos,
                   data_size);
        }
    }
}
