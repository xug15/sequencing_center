# More code for Biomarker liquid

## Create table browseebd

```sh
mysql -uroot -padmin
use BiomarkerLiquid;

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

```sh
CREATE TABLE biomarkersearch
(
Disease VARCHAR(50), 
Specimen VARCHAR(50), 
Datasets VARCHAR(50), 
Level VARCHAR(50), 
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
INDEX name (Disease,Specimen,Datasets)
) ENGINE=MyISAM ;
```

```sh
#mysqlimport -u root -padmin -d --local --fields-terminated-by="\t" --lines-terminated-by="\n" database table.sql;

mysqlimport -u root -padmin -d --local --fields-terminated-by="\t" --lines-terminated-by="\n" BiomarkerLiquid browseebd.sql;
mysqlimport -u root -padmin -d --local --fields-terminated-by="\t" --lines-terminated-by="\n" BiomarkerLiquid biomarkersearch.sql;

```

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
        $result = $db->query($query);
        $resultArray = $result->fetchAll();
        return $resultArray;
}
$response=gtfinfor($disease,$specimen,$datasets,$table);
echo json_encode($response);
?>

```

**browse.php**

+ http://lulab.life.tsinghua.edu.cn/ebd/php/browse.php?disease=Colorectal cancer&specimen=Plasma EV&datasets=PRJNA540919

```php
<?php
ini_set('memory_limit', '1024M');
$db = new PDO('mysql:host=127.0.0.1;dbname=BiomarkerLiquid', 'root', 'admin');
$table='browseebd';

$msgArray = array('code'=>0, 'data'=>array(), 'message'=>'参数接收错误，请关闭浏览器后重试。');
$disease = isset($_POST['disease']) ? trim($_POST['disease']) : trim($_GET['disease']);
$specimen = isset($_POST['specimen']) ? trim($_POST['specimen']) : trim($_GET['specimen']);
$datasets = isset($_POST['datasets']) ? trim($_POST['datasets']) : trim($_GET['datasets']);

header('content-type:application:json;charset=utf8');
header('Access-Control-Allow-Origin:*');
header('Access-Control-Allow-Methods:POST');
header('Access-Control-Allow-Headers:x-requested-with,content-type');


function gtfinfor($disease,$specimen,$datasets){
        global $db, $table;
        $query = "select * from ".$table." where Disease='".$disease."' and Specimen ='".$specimen."' and Datasets='".$datasets."'";
        $result = $db->query($query);
        $resultArray = $result->fetchAll();
        return $resultArray;
}
$response=gtfinfor($disease,$specimen,$datasets);
echo json_encode($response);
?>

```


