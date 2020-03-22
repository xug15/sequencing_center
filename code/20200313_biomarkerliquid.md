# MySQL

## create database.
```sh
CREATE DATABASE BiomarkerLiquid CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
use BiomarkerLiquid;
```
## create table
```sh
CREATE TABLE browseebd
(
Disease VARCHAR(50),
Specimen VARCHAR(50),
Datasets VARCHAR(50),
Genesymbol VARCHAR(50),
RNAtype VARCHAR(50),
Detectionfreq FLOAT,
Cancerexp FLOAT,
Normalexp FLOAT,
log2FC FLOAT,
FDR FLOAT,
ROCAUC FLOAT,
PanelAUC FLOAT,
Biomarkerreports VARCHAR(50),
sROCAUC FLOAT,
Markerpoteintial VARCHAR(50),
Summary TEXT,
INDEX name (Disease,Specimen,Datasets)
) ENGINE=MyISAM ;
```
|type |description |
|-| -|
|CHAR(M)| − A fixed-length string between 1 and 255 characters in length (for example CHAR(5)), right-padded with spaces to the specified length when stored. Defining a length is not required, but the default is 1. |
|VARCHAR(M)| − A variable-length string between 1 and 255 characters in length. For example, VARCHAR(25). You must define a length when creating a VARCHAR field. |
|TEXT |− A field with a maximum length of 65535 characters. BLOBs are "Binary Large Objects" and are used to store large amounts of binary data, such as images or other types of files. Fields defined as TEXT also hold large amounts of data. The difference between the two is that the sorts and comparisons on the stored data are case| 
|INT[Length]|Range of –2,147,483,648 to 2,147,483,647 or 0 to 4,294,967,295|
|BIGINT[Length]|Range of –9,223,372,036,854,775,808 to 9,223,372,036,854,775,807 or 0 to 18,446,744,073,709,551,615 unsigned|
|FLOAT[Length, Decimals]|A small number with a floating decimal point|

[more code](./)

```sh
show tables;
```

## import data
```sh
# Description:
mysqlimport -u root -padmin -d --local --fields-terminated-by="\t" --lines-terminated-by="\n" database table;


mysqlimport -u root -padmin -d --local --fields-terminated-by="\t" --lines-terminated-by="\n" BiomarkerLiquid browseebd;

```

## test mysql

```sh
select * from browseebd where Disease="Colorectal cancer" and Specimen ="Plasma EV" and Datasets="PRJNA540919";
select * from biomarkersearch where Disease="Colorectal cancer" and Specimen ="Plasma EV" and Datasets="PRJNA540919";

```

## php code. 

**biomarkersearch.php**

+ http://lulab.life.tsinghua.edu.cn/ebd/php/biomarkersearch.php?disease=Colorectal cancer&specimen=Plasma EV&datasets=PRJNA540919


```php
<?php
ini_set('memory_limit', '1024M');
$db = new PDO('mysql:host=127.0.0.1;dbname=BiomarkerLiquid', 'root', 'admin');
$table='biomarkersearch';

$msgArray = array('code'=>0, 'data'=>array(), 'message'=>'参数接收错误，请关闭浏览器后重试。');
$disease = isset($_POST['disease']) ? trim($_POST['disease']) : trim($_GET['disease']);
$specimen = isset($_POST['specimen']) ? trim($_POST['specimen']) : trim($_GET['specimen']);
$datasets = isset($_POST['datasets']) ? trim($_POST['datasets']) : trim($_GET['datasets']);

header('content-type:application:json;charset=utf8');
header('Access-Control-Allow-Origin:*');
header('Access-Control-Allow-Methods:POST');
header('Access-Control-Allow-Headers:x-requested-with,content-type');


function gtfinfor($disease,$specimen,$datasets,$table){
        global $db ;
        $query = "select * from ".$table." where Disease='".$disease."' and Specimen ='".$specimen."' and Datasets='".$datasets."'";
        //$query="select * from biomarkersearch where Disease='Colorectal cancer' and Specimen ='Plasma EV' and Datasets='PRJNA540919'";
        //$query="select * from gtf where transcript_index='ENST00000516445'";
        $result = $db->query($query);
        $resultArray = $result->fetchAll();
        return $resultArray;
}
$response=gtfinfor($disease,$specimen,$datasets,$table);
echo json_encode($response);
?>
```

## html

```html
    <style>
      .loader {
  border: 16px solid #f3f3f3;
  border-radius: 50%;
  border-top: 16px solid #3498db;
  border-bottom: 16px solid #3498db;
  width: 120px;
  height: 120px;
  -webkit-animation: spin 2s linear infinite;
  animation: spin 2s linear infinite;
}

@-webkit-keyframes spin {
  0% { -webkit-transform: rotate(0deg); }
  100% { -webkit-transform: rotate(360deg); }
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}
    </style>

    <!-- Insert this line after script imports -->
    <script type="text/javascript">
      if (window.module)
        module = window.module;
    </script>

<select id="sel_disease" style="font-size: 18px ;margin-top: 0%; float: inherit; width: 100%; height: 50px; text-align: center;">
              <option value="Colorectal cancer">Colorectal cancer</option>
              <option value="Liver cancer">Liver cancer</option>
              <option value="Pancreatic_cancer">Pancreatic cancer</option>
              
              <option value="Non_small_cell_lung_carcinoma">Non small cell lung carcinoma</option>
              <option value="Hepatic_cell_carcinoma">Hepatic cell carcinoma</option>
              <option value="Glioblastoma">Glioblastoma</option>
              <option value="Breast_cancer">Breast cancer</option>
              <option value="Serous">Serous</option>
              <option value="Chronic_fatigue_syndrome">Chronic fatigue syndrome</option>
            </select>

<div class="row" style="margin:10px 0; display:none" id="loader_box">
      <div style="width: 14%;margin: 0 auto;">
          <div class="loader"></div>
      </div>

<div id="box_result" class="row row_con" style="display:none;">
      <div class="col-md-12 col-xs-12">
          <div class="current_setting head_title" style="margin-bottom:10px;">
            <p id='title_result'>
              Current settings:
              <br />
              Disease: Liver cancer; Specimen: Plasma EV; Dataset(s): GSE100207
            </p>
          </div>
          <div class="row">
            <div id="information_result"></div>
            
          </div>
      </div>
    </div>
```


## javascript

```js
   infor_content='';
  function get_information(disease,specimen,datasets)
{
//var key="ENST00000384052";

var content=$.ajax({ 

    type:'GET',
    dataType:'json',
    url:'php/biomarkersearch.php',
    data:'disease='+disease+'&specimen='+specimen+'&datasets='+datasets,
    error:function(jqXHR, textStatus, errorThrown)
    {
        //alert("jqXHR:"+jqXHR+"\ntextStatus:"+textStatus+"\nerrorThrown:"+errorThrown);
        console.log("jqXHR:"+jqXHR+"\ntextStatus:"+textStatus+"\nerrorThrown:"+errorThrown);
    },
    success:function(response)
        {
          data_response=response;
          if(response.length!==0){
            infor_content+='<table id="example2" class="display" style="width:100%"><thead><tr><th>Gene symbol	</th><th>RNA type</th><th>Detection freq.	</th><th>Cancer exp.</th><th>Normal exp.</th><th>logFC</th>';
            infor_content+='<th>FDR</th><th>ROC AUC</th><th>Panel AUC</th><th>Biomarker reports</th><th>Meta-analysis (sROC AUC)</th><th>Summary</th></tr></thead><body>'
            for (var i=0;i<response.length;i++)
            {
              infor_content+='<tr><th>'+data_response[i][4]+'</th><th>'+data_response[i][5]+'</th><th>'+data_response[i][6]+'</th><th>'+data_response[i][7]+'</th><th>'+data_response[i][8]+'</th><th>'+data_response[i][9]+'</th><th>'+data_response[i][10]+'</th><th>'+data_response[i][11]+'</th><th>'+data_response[i][12]+'</th><th>'+data_response[i][13]+'</th><th>'+data_response[i][14]+'</th><th>'+data_response[i][15]+'</th></tr>';
            }
            infor_content+='</body></table>';
          document.getElementById('information_result').innerHTML = infor_content;
          $('#example2').DataTable();
          document.getElementById("loader_box").style.display = "none";
          document.getElementById("box_result").style.display = "block";
          
          }
          
          
        //console.log(response);
        //alert(infor_content);
        //alert(infor_content2);
        return response;
    }
    })
    return content;
}

    function  search()
    {
      
      document.getElementById("loader_box").style.display = "block";
      disease=document.getElementById('sel_disease').value;
      specimen=document.getElementById('specimen').value;
      datasets=document.getElementById('datasets').value;
      get_information(disease,specimen,datasets);
      document.getElementById("box_input").style.display = "none";
      
      document.getElementById('title_result').innerHTML = 'Current settings:<br />Disease: '+disease+'; Specimen: '+specimen+'; Dataset(s): '+datasets;
      console.log(disease);
      console.log(specimen);
      console.log(datasets);
    }

```


