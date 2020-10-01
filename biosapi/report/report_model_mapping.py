from biosapi import BiosModelMapper

def report_species_mapping(target_database, min_score, model_ids, api):
    report = {
        'models' : list(model_ids),
        'databases' : [target_database],
        'records' : [],
        'data' : {}
    }
    all_data = {}
    database_id_to_model = {}
    for model_id in model_ids:
        all_data[model_id] = {}
        model_mapper = BiosModelMapper(api, model_id)
        mapping = model_mapper.get_spi_annotation(target_database, min_score)
        compartments = {}
        for o in model_mapper.cmp:
            compartments[o['id']] = o['bios_scmp_entry'] if 'bios_scmp_entry' in o else f"{o['id']}@{model_id}"
        for spi_id in mapping:
            uid = model_mapper.spi_sid_to_uid[spi_id]
            spi = model_mapper.spi[uid]
            all_data[model_id][spi_id] = {
                'bios_id' : spi['bios_id'],
                'name' : spi['name']
            }
            database_id = mapping[spi_id]
            if database_id not in database_id_to_model:
                database_id_to_model[database_id] = {}
            spi_cmp = compartments[spi['compartment']]
            if spi_cmp not in database_id_to_model[database_id]:
                database_id_to_model[database_id][spi_cmp] = {}
            if model_id not in database_id_to_model[database_id][spi_cmp]:
                database_id_to_model[database_id][spi_cmp][model_id] = []
            database_id_to_model[database_id][spi_cmp][model_id].append(spi_id)
    report['data'] = all_data
    for cpd_id in database_id_to_model:
        record = {
            'model_data' : database_id_to_model[cpd_id],
            'database': {target_database: [cpd_id]}
        }
        report['records'].append(record)
    
    return report