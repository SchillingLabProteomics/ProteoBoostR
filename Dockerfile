# escape=`

ARG BASE_IMAGE=mcr.microsoft.com/windows:ltsc2019
FROM ${BASE_IMAGE}

SHELL ["powershell", "-Command", "$ErrorActionPreference = 'Stop';"]
RUN Invoke-WebRequest -Uri "https://cran.r-project.org/bin/windows/base/old/4.3.3/R-4.3.3-win.exe" -OutFile "R-4.3.3-win.exe"; `
    Start-Process -Wait -FilePath "R-4.3.3-win.exe" -ArgumentList '/VERYSILENT','/SUPPRESSMSGBOXES','/NORESTART','/SP-','/DIR=C:\R\R-4.3.3'; `
    Remove-Item R-4.3.3-win.exe

RUN Invoke-WebRequest -Uri "https://cran.r-project.org/bin/windows/Rtools/rtools43/files/rtools43-5976-5975.exe" -OutFile "Rtools43.exe"; `
    Start-Process -Wait -FilePath "Rtools43.exe" -ArgumentList '/VERYSILENT','/SUPPRESSMSGBOXES','/NORESTART','/SP-','/DIR=C:\Rtools43'; `
    Remove-Item Rtools43.exe

COPY app.R C:\shinyapp\app.R
COPY install_requirements.R C:\shinyapp\install_requirements.R
RUN & "'C:\R\R-4.3.3\bin\Rscript.exe'" "C:\shinyapp\install_requirements.R"
RUN Write-Output 'C:\shinyapp\outputs' | Out-File -FilePath C:\shinyapp\ProteoBoostR.filesystem -Encoding ASCII
RUN New-Item -ItemType Directory -Path 'C:\shinyapp\outputs' -Force
EXPOSE 3838
CMD ["C:\\R\\R-4.3.3\\bin\\Rscript.exe", "-e", "shiny::runApp('C:/shinyapp/app.R', host='0.0.0.0', port=3838)"]