drop table if exists tblesv;
create table tblesv (
  centroid_id integer primary key,
  centroid text not null,
  lat text default '',
  lon text default '',
  depth text default '', 
  relative_abundance text default '',
  esv_tempreature text default '',
  esv_salinity text default '',
  cruise_name text default '',
  size_frac_lower real default 0.0,
  size_frac_upper real default 0.0,
  id text default ''
);

create index centroid on tblesv (centroid);
